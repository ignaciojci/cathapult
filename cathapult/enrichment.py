import pandas as pd
import numpy as np
from .analyze import extract_domain_levels, annotate_domains, deep_annotate_domains
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def collapse_domain_levels(df):
    """
    Optimized: Uses vectorized melting instead of manual iteration.
    """
    df = df[df["cath_label"] != "-"].copy()
    df["domain"] = df["cath_label"]
    df = extract_domain_levels(df)
    
    domain_cols = ["domain.first.level", "domain.two.levels", "domain.three.levels", "domain"]
    
    # Vectorized reshape (melt) is significantly faster than iterrows
    df_features = df.melt(
        id_vars=['uniprot_acc'], 
        value_vars=domain_cols, 
        var_name='level_name', 
        value_name='domain'
    )
    
    # Map level names back to numeric indices for compatibility
    level_map = {col: i for i, col in enumerate(domain_cols)}
    df_features['level'] = df_features['level_name'].map(level_map)
    
    return df_features[['domain', 'level', 'uniprot_acc']]

def calculate_odds_ratio(df1, df2, unique_features=False, min_count=5):
    """
    Optimized: Uses value counts pre-calculation and threshold filtering.
    """
    df1_features = collapse_domain_levels(df1)
    df2_features = collapse_domain_levels(df2)
    
    if unique_features:
        df1_features = df1_features.drop_duplicates(subset=['domain','uniprot_acc'])
        df2_features = df2_features.drop_duplicates(subset=['domain','uniprot_acc'])
    
    
    # Combine unique features
    all_results = pd.DataFrame()
    
    for level in [0, 1, 2, 3]:
        results = []
        df1_features_level = df1_features[df1_features['level'] == level]
        df2_features_level = df2_features[df2_features['level'] == level]

        # Pre-calculate counts once
        counts1 = df1_features_level['domain'].value_counts()
        counts2 = df2_features_level['domain'].value_counts()

        all_features = set(counts1.index) | set(counts2.index)
        
        total_group1 = len(df1_features_level)
        total_group2 = len(df2_features_level)

        for feature in all_features:
            a = counts1.get(feature, 0)
            c = counts2.get(feature, 0)
            
            b = total_group1 - a
            d = total_group2 - c

            # --- NEW: Skip rare levels ---
            if a < min_count or c < 1:
                results.append({
                    'feature': feature,
                    'grp1_count': a, 'grp1_total': total_group1,
                    'grp2_count': c, 'grp2_total': total_group2,
                    'odds.ratio': np.nan,
                    'log.odds.ratio': np.nan,
                    'p.value': np.nan,
                    'ci.lower': np.nan,
                    'ci.upper': np.nan,
                    'domain': feature,
                    'domain.level': level+1,
                })
                continue
                
            # Contingency table calculation
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
            
            # Handling Confidence Intervals
            if a == 0 or b == 0 or c == 0 or d == 0:
                ci_lower, ci_upper = np.nan, np.nan
            else:
                # Woolf's method
                log_or = np.log(odds_ratio)
                se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
                ci_lower = np.exp(log_or - 1.96 * se_log_or)
                ci_upper = np.exp(log_or + 1.96 * se_log_or)

            results.append({
                'feature': feature,
                'grp1_count': a, 'grp1_total': total_group1,
                'grp2_count': c, 'grp2_total': total_group2,
                'odds.ratio': odds_ratio,
                'log.odds.ratio': np.log2(odds_ratio) if odds_ratio > 0 else -np.inf,
                'p.value': p_value,
                'ci.lower': np.log2(ci_lower) if not np.isnan(ci_lower) else np.nan,
                'ci.upper': np.log2(ci_upper) if not np.isnan(ci_upper) else np.nan,
                'domain': feature,
                'domain.level': level+1,
            })

        results_df = pd.DataFrame(results)
        
        if not results_df.empty:
            p_values = results_df['p.value']
            valid_idx = p_values.notna()

            # Initialize p.adj with NaN
            results_df['p.adj'] = np.nan

            if valid_idx.sum() > 0:
                _, p_adj, _, _ = multipletests(
                    p_values[valid_idx],
                    method='fdr_bh'
                )
                results_df.loc[valid_idx, 'p.adj'] = p_adj

            results_df = deep_annotate_domains(results_df)
            all_results = all_results.append(results_df)
            
    return all_results

def plot_odds_ratio(results_df, output_plot, alpha=0.05):
    """
    Generates and saves a plot of the odds ratios.
    """
    results_df = results_df.sort_values(by='log.odds.ratio', ascending=True).dropna(subset=['log.odds.ratio'])
    results_df['significant'] = results_df['p.value'] < alpha
    results_df['y.labels'] = results_df['feature'] + ' ' + results_df['domain.name']
    results_df = results_df.dropna()
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, max(6, len(results_df) * 0.275)))
    
    colors = {True: '#d62728', False: 'grey'} # Red for significant, grey for not
    
    # Plot confidence intervals
    for _, row in results_df.iterrows():
        ax.plot([row['ci.lower'], row['ci.upper']], [row['y.labels'], row['y.labels']], 'k-', lw=1, zorder=1)
    

    # Plot odds ratio points
    ax.scatter(results_df['log.odds.ratio'], results_df['y.labels'], 
               c=results_df['significant'].map(colors),
               s=60, zorder=2, edgecolors='black', linewidths=0.7)

    ax.axvline(x=1, color='black', linestyle='--', lw=1)
    #ax.set_xscale('log2')
    ax.set_xlabel('Odds Ratio (log2 scale)\n(Group 1 vs Group 2)', fontsize=12)
    ax.set_ylabel('Feature', fontsize=12)
    ax.set_title('Feature Enrichment Odds Ratio', fontsize=14, weight='bold')
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.3, right=0.98, top=0.95, bottom=0.15)
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    logging.info(f"📊 Plot saved to {output_plot}")