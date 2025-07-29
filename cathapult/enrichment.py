import pandas as pd
import numpy as np
from .analyze import extract_domain_levels, annotate_domains
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def collapse_domain_levels(df):
    df = df[df["cath_label"] != "-"]
    df["domain"] = df["cath_label"]
    df = extract_domain_levels(df)
    
    domain_cols = ["domain.first.level", "domain.two.levels", "domain.three.levels", "domain"]
    
    df_features = pd.DataFrame(columns=['domain', 'level', 'uniprot_acc'])
    iout = 0
    
    for irow, row in df.iterrows():
        for icol, col in enumerate(domain_cols):
            df_features.loc[iout] = [row[col], icol, row['uniprot_acc']]
            iout = iout + 1
            
    return(df_features)
    

def calculate_odds_ratio(df1, df2, unique_features=False):
    """
    Calculates the odds ratio for each feature between two groups.
    """
    
    df1_features = collapse_domain_levels(df1)
    df2_features = collapse_domain_levels(df2)
    
    if unique_features:
        df1_features = df1_features.drop_duplicates(subset=['domain','uniprot_acc'])
        df2_features = df2_features.drop_duplicates(subset=['domain','uniprot_acc'])
    
    all_features = sorted(list(set(df1_features['domain'].unique()) | set(df2_features['domain'].unique())))
    results = []

    total_group1 = len(df1)
    total_group2 = len(df2)
    
    for feature in all_features:
        # Contingency table: [[group1_has, group1_not_has], [group2_has, group2_not_has]]
        a = df1_features['domain'].value_counts().get(feature, 0)
        b = total_group1 - a
        c = df2_features['domain'].value_counts().get(feature, 0)
        d = total_group2 - c
        
        # Fisher's exact test requires non-zero entries for CI calculation
        # We handle zero counts gracefully
        if a == 0 or b == 0 or c == 0 or d == 0:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
            ci_lower, ci_upper = np.nan, np.nan
        else:
            table = np.array([[a, b], [c, d]])
            odds_ratio, p_value = fisher_exact(table)
            # Woolf's method for confidence interval of log(OR)
            log_or = np.log(odds_ratio)
            se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)

        results.append({
            'feature': feature,
            'grp1_count' : a,
            'grp1_total': b,
            'grp2_count' : c,
            'grp2_total' : d,
            'grp1_proportion' : a/b,
            'grp2_proportion' : c/d,
            'odds.ratio' : odds_ratio,
            'log.odds.ratio': np.log2(odds_ratio),
            'p.value': p_value,
            'ci.lower': np.log2(ci_lower),
            'ci.upper': np.log2(ci_upper),
            'domain': feature,
        })

    results_df = pd.DataFrame(results)

    if not results_df.empty and 'p.value' in results_df:
        # Benjamini-Hochberg p-value correction
        _, p_adj, _, _ = multipletests(results_df['p.value'].dropna(), alpha=0.05, method='fdr_bh')
        results_df.loc[results_df['p.value'].notna(), 'p.adj'] = p_adj

    results_df = annotate_domains(results_df)    
    
    return results_df

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