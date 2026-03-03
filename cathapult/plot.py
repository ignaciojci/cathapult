import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(
    "/broad/grekalab/Peach Arines/Data/Proteomics/HuProt_TMED9-GOLD-CC/2025 January Re-run/Analyses/TED Domains/secretome_odds_ratio/odds_ratio_6k_hp_log.tsv",
    sep="\t"
)

# -----------------------------
# Column headers (configurable)
# -----------------------------
p_value_header = "p.value"
odds_ratio_header = "odds.ratio"
feature_id_header = "domain"
feature_name_header = "domain.name"
count_header = "grp1_count"
group_header = "domain.level"
feature_name_headers = [
    "domain.name",
    "domain.third.level.name",
    "domain.second.level.name",
    "domain.first.level.name",
]

# -----------------------------
# Parameters
# -----------------------------
TOP_N = 20
SIZE_SCALE = 0.4
LABEL_OFFSET_FRAC = 0.08
Y_SPACING = 0.8
ROW_HEIGHT = 0.45
MIN_HEIGHT = 2.5
MAX_HEIGHT = 8

# -----------------------------
# Prepare data (using headers)
# -----------------------------
df["resolved_feature_name"] = (
    df[feature_name_headers]
    .bfill(axis=1)
    .iloc[:, 0]
)

for group_value in [2, 3, 4]:
    df_filtered = df[df[group_header] == group_value]
    df_filtered = df_filtered[df[odds_ratio_header] > 1]
    
    TOP_N_group = min(TOP_N, len(df))
    
    df_plot = (
        df_filtered
        .copy()
        .sort_values(p_value_header)
        .head(TOP_N_group)
    )
    
    df_plot["neg_log10_p"] = -np.log10(df_plot[p_value_header])
    
    # Reverse order so best hit is on top
    df_plot = df_plot.iloc[::-1].reset_index(drop=True)
    
    COUNT_MIN_SIZE = 30
    COUNT_MAX_SIZE = 300
    
    df_plot["log_count"] = np.log10(df_plot[count_header].clip(lower=1))
    
    # Normalize to [0, 1]
    log_count_norm = (
        (df_plot["log_count"] - df_plot["log_count"].min()) /
        (df_plot["log_count"].max() - df_plot["log_count"].min())
    )
    
    df_plot["point_size"] = (
        COUNT_MIN_SIZE +
        log_count_norm * (COUNT_MAX_SIZE - COUNT_MIN_SIZE)
    )
    
    y_pos = np.arange(len(df_plot)) * Y_SPACING
    
    # -----------------------------
    # Plot
    # -----------------------------
    fig_height = min(
        MAX_HEIGHT,
        max(MIN_HEIGHT, ROW_HEIGHT * TOP_N_group + 1)
    )
    fig, ax = plt.subplots(figsize=(7, fig_height))
    
    sc = ax.scatter(
        df_plot["neg_log10_p"],
        y_pos,
        s=df_plot["point_size"],
        c=df_plot[odds_ratio_header],
        cmap="viridis",
        norm=LogNorm(
            vmin=df_plot[odds_ratio_header].min(),
            vmax=df_plot[odds_ratio_header].max()
        ),
        edgecolor="black",
        linewidth=0.5,
        zorder=3
    )
    
    # -----------------------------
    # Axes formatting
    # -----------------------------
    ax.set_yticks(y_pos)
    ax.set_yticklabels([])
    ax.set_xlabel(r"$- \log_{10}(\mathrm{p\ value})$", fontsize=12)
    ax.set_title("Top enriched features", fontsize=14, weight="bold")
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=11)
    ax.tick_params(axis="y", length=0)
    
    # -----------------------------
    # Feature labels (right of dots)
    # -----------------------------
    x_min, x_max = ax.get_xlim()
    x_offset = (x_max - x_min) * LABEL_OFFSET_FRAC
    
    for i, row in df_plot.iterrows():
        # Bold feature ID
        ax.text(
            row["neg_log10_p"] + x_offset,
            y_pos[i],
            str(row[feature_id_header]),
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold" 
        )
        
        # Feature name (normal weight)
        ax.text(
            row["neg_log10_p"] + x_offset,
            y_pos[i] - 0.4 * Y_SPACING,
            str(row["resolved_feature_name"]),
            va="center",
            ha="left",
            fontsize=9
        )
        
    # Expand x-limits so labels fit
    ax.set_xlim(x_min, x_max + (x_max - x_min) * 2.5)
    #ax.set_ylim(1, 7)
    
    # -----------------------------
    # Colorbar (odds ratio)
    cbar = plt.colorbar(
        sc,
        ax=ax,
        pad=0.02,
        shrink=0.3,      # <-- makes the bar shorter
        aspect=20,        # <-- keeps it slim and Prism-like
        format="%d"
    )
    
    cbar.set_label("Odds ratio", fontsize=11)
    cbar.ax.tick_params(labelsize=10)
    
    # -----------------------------
    # Size legend (count)
    # -----------------------------
    legend_counts = np.percentile(
        df_plot[count_header],
        [25, 50, 75]
    ).astype(int)
    
    for c in legend_counts:
        size = (
            COUNT_MIN_SIZE +
            (np.log10(c) - df_plot["log_count"].min()) /
            (df_plot["log_count"].max() - df_plot["log_count"].min())
            * (COUNT_MAX_SIZE - COUNT_MIN_SIZE)
        )
        ax.scatter(
            [],
            [],
            s=size,
            edgecolors="black",
            facecolors="none",
            label=str(c)
        )
        
    ax.legend(
        title="Observation count",
        frameon=False,
        loc="lower right",
        fontsize=10,
        title_fontsize=11
    )
    
    plt.tight_layout()
    plt.show()
    plt.savefig(f"plot_group_" + str(group_value) + ".png")

