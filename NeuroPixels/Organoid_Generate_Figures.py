#!/usr/bin/env python3
"""
Organoid 1/f aperiodic exponent: cNS vs cNS+GBM.
Loads per-channel FOOOF slope CSVs, aggregates to sample level,
and plots Tumor co-culture vs Neuron Only boxplot with stats.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_CSV = os.path.join(BASE_DIR, "organoid_csvs", "organoid_oneoverf_all_channels.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "organoid_figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROUP_COLORS = {
    "neuron_only": "#3274A1",
    "tumor":       "#E1812C",
}


def plot_group_boxplot(df, agg_col, level_label):
    """Boxplot: Tumor co-culture vs Neuron Only with rank-sum p-value."""
    agg = df.groupby([agg_col, 'Group'])['Slope'].mean().reset_index()

    neuron_vals = agg.loc[agg['Group'] == 'neuron_only', 'Slope'].values
    tumor_vals = agg.loc[agg['Group'] == 'tumor', 'Slope'].values

    order = ["neuron_only", "tumor"]
    data_dict = {"neuron_only": neuron_vals, "tumor": tumor_vals}
    xlabel_map = {"neuron_only": "cNS", "tumor": "cNS + GBM"}

    fig, ax = plt.subplots(figsize=(5, 6))
    positions = [1, 2]

    for i, label in enumerate(order):
        vals = data_dict[label]
        if len(vals) == 0:
            continue
        pos = positions[i]
        bp = ax.boxplot(
            [vals], positions=[pos], widths=0.6,
            patch_artist=True, showfliers=False,
        )
        bp['boxes'][0].set_facecolor(GROUP_COLORS[label])
        bp['boxes'][0].set_alpha(0.6)
        jitter = np.random.uniform(-0.12, 0.12, len(vals))
        ax.scatter(
            np.full(len(vals), pos) + jitter, vals,
            color=GROUP_COLORS[label], edgecolor='black',
            s=50, alpha=0.8, zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels([xlabel_map[l] for l in order], fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    if len(neuron_vals) > 0 and len(tumor_vals) > 0:
        _, p = ranksums(neuron_vals, tumor_vals)
        y_max = max(neuron_vals.max(), tumor_vals.max())
        y_range = y_max - min(neuron_vals.min(), tumor_vals.min())
        y_bar = y_max + 0.08 * y_range
        ax.plot([1, 2], [y_bar, y_bar], 'k-', linewidth=0.8)
        star = '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
        ax.text(1.5, y_bar + 0.02 * y_range, star,
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_ylabel('Aperiodic Exponent', fontsize=12)
    ax.set_title(f'Neuropixels Organoid: 70-150 Hz\n({level_label})',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    fname = f'Organoid_OneOverF_cNS_vs_GBM_{level_label.replace(" ", "_")}.pdf'
    out_path = os.path.join(OUTPUT_DIR, fname)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    print("=" * 60)
    print("ORGANOID 1/f: cNS vs cNS+GBM")
    print("=" * 60)

    if not os.path.isfile(DATA_CSV):
        print(f"ERROR: Combined CSV not found at {DATA_CSV}")
        print("Run Organoid_Compute_FOOOF.py first.")
        exit(1)

    df = pd.read_csv(DATA_CSV)

    gbm_samples = [s for s, m in {
        "HCN-1": "neuron_only", "HCN-2": "neuron_only", "HCN-3": "neuron_only",
        "SF725-1": "tumor", "SF725-2": "tumor", "SF725-3": "tumor", "SF725-4": "tumor",
    }.items()]
    df = df[df['Sample'].isin(gbm_samples)]
    print(f"Loaded {len(df)} channel-level rows (GBM + Neuron Only)")

    plot_group_boxplot(df, agg_col='Sample', level_label='Sample Level')

    print("\nDone!")
