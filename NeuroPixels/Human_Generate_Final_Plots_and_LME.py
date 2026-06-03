#!/usr/bin/env python3
"""
Aperiodic exponent by glioma subtype (Oligo, Astro, GBM) across cortical depth.
Bar+Scatter plots: empty bars with outline, dots on top, mean +/- SEM.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, ranksums
import statsmodels.stats.multitest as smm

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

data_dir = os.path.join(BASE_DIR, 'patient_csvs')
output_dir = os.path.join(BASE_DIR, 'output', 'figures')
os.makedirs(output_dir, exist_ok=True)

pial_surfaces = {
    "NP87_B1": 6820,
    "NP32_B2": 5930, "NP32_B3": 5910, "NP37_B1": 6350, "NP37_B2": 6350,
    "NP38_B5": 6430, "NP38_B6": 5830, "NP46_B1": 6180, "NP46_B2": 6917,
    "NP54_B1": 6670, "NP64_B1": 7440, "NP65_B2": 7500, "NP78_B1": 7550,
    "NP93_B1": 5000, "NP94_B1": 6639, "NP94_B3": 6000, "NP95_B1": 6600,
    "NP119_B1": 7000, "NP124_B3": 5658, "NP125_B1": 6473,
    "NP123_B1": 7500, "NP130_B2": 6637
}

tumor_patients = ["NP32", "NP37", "NP38", "NP46", "NP54", "NP64", "NP65",
                  "NP78", "NP93", "NP94", "NP95", "NP123", "NP130"]

subtype_map = {
    "NP32": "Oligo", "NP38": "Oligo",
    "NP37": "Astro", "NP78": "Astro", "NP93": "Astro", "NP95": "Astro", "NP123": "Astro",
    "NP46": "GBM", "NP54": "GBM", "NP64": "GBM", "NP65": "GBM", "NP94": "GBM", "NP130": "GBM"
}

colors_subtype = {
    'Oligo': '#000000',
    'Astro': '#E91E8C',
    'GBM': '#00A99D'
}

depth_step = 1000
depth_bins = np.arange(0, 5001, depth_step)
depth_labels = [f'{d}-{d+depth_step}' for d in depth_bins[:-1]]


def load_tumor_data():
    all_data = []
    for patient in tumor_patients:
        patient_files = [f for f in os.listdir(data_dir)
                        if f.startswith(f"{patient}_") and f.endswith("_fooof_slope_by_channel.csv")]
        for file_name in patient_files:
            try:
                file_data = pd.read_csv(os.path.join(data_dir, file_name))
                parts = file_name.split('_')
                if len(parts) >= 3 and parts[1].startswith('B'):
                    block = parts[1]
                    patient_block_key = f"{patient}_{block}"
                    imec_id = file_name.split('imec')[1].split('_')[0] if 'imec' in file_name else '0'
                    probe_id = f"{patient}_{block}_imec{imec_id}"
                else:
                    continue
                if patient_block_key not in pial_surfaces:
                    continue
                pial_value = pial_surfaces[patient_block_key]
                file_data["Distance_from_surface"] = pial_value - file_data["Distance_to_Tip"]
                file_data = file_data[(file_data['Distance_from_surface'] >= 0) & (file_data['Distance_from_surface'] <= 5000)]
                if len(file_data) == 0:
                    continue
                file_data['Depth_Bin'] = pd.cut(file_data['Distance_from_surface'], bins=depth_bins, labels=depth_labels, right=False)
                file_data['Patient_ID'] = patient
                file_data['Probe_ID'] = probe_id
                file_data['Subtype'] = subtype_map.get(patient, 'Unknown')
                all_data.append(file_data)
            except Exception:
                continue
    return pd.concat(all_data, ignore_index=True).dropna(subset=['Slope', 'Depth_Bin'])


def create_subtype_bar_scatter(df, level_name):
    subtype_order = ['Oligo', 'Astro', 'GBM']
    labels = ['O', 'A', 'GBM']

    n_bins = len(depth_labels)
    fig, axes = plt.subplots(1, n_bins, figsize=(3.0 * n_bins, 5), sharey=True)

    for i, depth_bin in enumerate(depth_labels):
        ax = axes[i]
        positions = list(range(1, len(subtype_order) + 1))
        bar_width = 0.6

        for j, subtype in enumerate(subtype_order):
            data = df[(df['Subtype'] == subtype) & (df['Depth_Bin'] == depth_bin)]['Slope'].dropna().values

            if len(data) > 0:
                mean_val = np.mean(data)
                sem_val = sem(data) if len(data) > 1 else 0
                color = colors_subtype[subtype]

                ax.bar(positions[j], mean_val, width=bar_width,
                      facecolor='white', edgecolor=color, linewidth=2)

                ax.errorbar(positions[j], mean_val, yerr=sem_val,
                           fmt='none', ecolor='black', capsize=3, capthick=1.5, elinewidth=1.5)

                jitter = np.random.uniform(-0.2, 0.2, len(data))
                ax.scatter([positions[j]] * len(data) + jitter, data,
                          color=color, edgecolor=color, linewidth=0.5,
                          s=15, alpha=0.7, zorder=3)

        # Pairwise rank-sum tests with FDR correction
        bin_data = {s: df[(df['Subtype'] == s) & (df['Depth_Bin'] == depth_bin)]['Slope'].dropna().values
                    for s in subtype_order}
        pairs = [(0, 1), (0, 2), (1, 2)]
        raw_pvals = []
        valid_pairs = []
        for a, b in pairs:
            d1, d2 = bin_data[subtype_order[a]], bin_data[subtype_order[b]]
            if len(d1) > 0 and len(d2) > 0:
                _, p = ranksums(d1, d2)
                raw_pvals.append(p)
                valid_pairs.append((a, b))

        if raw_pvals:
            _, fdr_pvals, _, _ = smm.multipletests(raw_pvals, method='fdr_bh')
            all_vals = np.concatenate([v for v in bin_data.values() if len(v) > 0])
            y_max = all_vals.max()
            y_range = y_max - all_vals.min() if len(all_vals) > 1 else 1
            bracket_idx = 0
            for (a, b), fdr_p in zip(valid_pairs, fdr_pvals):
                if fdr_p < 0.05:
                    star = '***' if fdr_p < 0.001 else '**' if fdr_p < 0.01 else '*'
                    y_bar = y_max + 0.08 * y_range * (bracket_idx + 1)
                    ax.plot([positions[a], positions[b]], [y_bar, y_bar],
                            'k-', linewidth=0.8)
                    ax.text((positions[a] + positions[b]) / 2, y_bar + 0.01 * y_range,
                            star, ha='center', va='bottom', fontsize=9, fontweight='bold')
                    bracket_idx += 1

        depth_mm = f"{int(depth_bin.split('-')[0])/1000:.0f}-{int(depth_bin.split('-')[1])/1000:.0f}mm"
        ax.set_title(depth_mm, fontsize=12, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=11)
        ax.set_xlim(0.3, len(subtype_order) + 0.7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel('Aperiodic Exponent', fontsize=12)

    fig.suptitle(f'Aperiodic Exponent by Glioma Subtype ({level_name} Level)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_file = os.path.join(output_dir, f'BarScatter_subtypes_{level_name.lower()}_level.pdf')
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_file}")


if __name__ == "__main__":
    print("Loading tumor patient data...")
    df_channel = load_tumor_data()

    df_probe = df_channel.groupby(['Patient_ID', 'Probe_ID', 'Subtype', 'Depth_Bin'],
                                   observed=True)['Slope'].mean().reset_index()
    df_patient = df_channel.groupby(['Patient_ID', 'Subtype', 'Depth_Bin'],
                                     observed=True)['Slope'].mean().reset_index()

    print(f"Data loaded: Channel={len(df_channel)}, Probe={len(df_probe)}, Patient={len(df_patient)}")

    print("\nGenerating subtype plots...")
    for level_name, df in [('Channel', df_channel), ('Probe', df_probe), ('Patient', df_patient)]:
        create_subtype_bar_scatter(df, level_name)

    print("\nDone!")
