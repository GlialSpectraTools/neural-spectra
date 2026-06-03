#!/usr/bin/env python3
"""
Organoid PSD overlay: GBM vs Neuron Only (mean +/- 95% CI).
Loads NWB files, computes channel-mean PSDs, caches as .npz,
then plots GBM and Neuron Only on a single axes.
"""

import os
import glob
import numpy as np
import pynapple as nap
from pynwb import NWBHDF5IO
from scipy.signal import welch, resample_poly
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = [
    'Nimbus Sans', 'Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 12

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NWB_BASE_PATH = os.path.join(BASE_DIR, "raw_nwb")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "organoid_figures")
CACHE_DIR = os.path.join(OUTPUT_DIR, "psd_cache")
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_FS = 400
LFP_DURATION_SEC = 180

SAMPLE_METADATA = {
    "HCN-1":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "HCN-2":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "HCN-3":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "SF725-1":  {"group": "tumor", "pathology": "Glioblastoma"},
    "SF725-2":  {"group": "tumor", "pathology": "Glioblastoma"},
    "SF725-3":  {"group": "tumor", "pathology": "Glioblastoma"},
    "SF725-4":  {"group": "tumor", "pathology": "Glioblastoma"},
}

PATHOLOGY_COLORS = {
    "Neuron Only":  "#BFBFBF",
    "Glioblastoma": "#7DBBB5",
}


def interpolate_notch_gaps(freqs, psd, notches=(60, 120, 180), tol=1.0):
    """Linearly interpolate across notch-filter dips in the PSD."""
    clean = psd.copy()
    for nf in notches:
        idx = np.where((freqs >= nf - tol) & (freqs <= nf + tol))[0]
        if len(idx) > 0:
            lo = max(0, idx[0] - 1)
            hi = min(len(freqs) - 1, idx[-1] + 1)
            clean[idx] = np.interp(
                freqs[idx], [freqs[lo], freqs[hi]],
                [clean[lo], clean[hi]])
    return clean


def discover_recordings():
    recordings = []
    for sample, meta in SAMPLE_METADATA.items():
        pattern = os.path.join(NWB_BASE_PATH, f"NPorganoid_{sample}_B*")
        dirs = sorted(glob.glob(pattern))
        for rec_dir in dirs:
            dirname = os.path.basename(rec_dir)
            if "plx" in dirname:
                continue
            nwb_file = os.path.join(rec_dir, f"{dirname}.nwb")
            if os.path.isfile(nwb_file):
                block = dirname.split("_B")[-1]
                recordings.append({
                    "sample": sample,
                    "block": f"B{block}",
                    "recording": dirname,
                    "nwb_path": nwb_file,
                    "group": meta["group"],
                    "pathology": meta["pathology"],
                })
    return recordings


def compute_psd(nwb_path, target_fs=TARGET_FS, duration_sec=LFP_DURATION_SEC):
    io = NWBHDF5IO(nwb_path, 'r')
    nwbfile = io.read()
    conversion = nwbfile.acquisition['ElectricalSeriesLF'].conversion
    io.close()

    data = nap.load_file(nwb_path)
    lfp_raw = data["ElectricalSeriesLF"]
    native_fs = lfp_raw.rate

    end_time = min(lfp_raw.time_support.start[0] + duration_sec,
                   lfp_raw.time_support.end[0])
    lfp_r = lfp_raw.restrict(
        nap.IntervalSet(start=lfp_raw.time_support.start[0], end=end_time)
    )
    lfp_np = np.asarray(lfp_r.values, dtype=np.float64) * conversion * 1e6

    lfp_ds = resample_poly(lfp_np, 4, 25, axis=0)

    nperseg = min(target_fs * 3, lfp_ds.shape[0])
    f, Pxx = welch(lfp_ds, fs=target_fs, nperseg=nperseg,
                   noverlap=nperseg // 2, axis=0)

    mean_psd = np.mean(Pxx, axis=1)
    return f, mean_psd


def load_all_psds(recordings):
    psd_data = {}
    for i, rec in enumerate(recordings):
        cache_path = os.path.join(CACHE_DIR, f"{rec['recording']}.npz")
        if os.path.isfile(cache_path):
            cached = np.load(cache_path)
            freqs, mean_psd = cached["freqs"], cached["mean_psd"]
            print(f"  [{i+1}/{len(recordings)}] {rec['recording']} (cached)")
        else:
            print(f"  [{i+1}/{len(recordings)}] {rec['recording']}...",
                  end=" ", flush=True)
            try:
                freqs, mean_psd = compute_psd(rec["nwb_path"])
                np.savez_compressed(cache_path, freqs=freqs, mean_psd=mean_psd)
                print("OK")
            except Exception as e:
                print(f"FAILED: {e}")
                continue
        mean_psd = interpolate_notch_gaps(freqs, mean_psd)
        psd_data[rec["recording"]] = (freqs, mean_psd)
    return psd_data


def plot_gbm_vs_neuron(recordings, psd_data):
    subset = ["Glioblastoma", "Neuron Only"]
    fig, ax = plt.subplots(figsize=(6, 5))
    handles = []
    for pathology in subset:
        color = PATHOLOGY_COLORS[pathology]
        recs = [r for r in recordings if r["pathology"] == pathology]
        psds = [psd_data[r["recording"]][1] for r in recs
                if r["recording"] in psd_data]
        if not psds:
            continue
        f = psd_data[recs[0]["recording"]][0]
        psd_stack = np.array(psds)
        mean = np.mean(psd_stack, axis=0)
        sem = np.std(psd_stack, axis=0) / np.sqrt(psd_stack.shape[0])
        ci = 1.96 * sem
        log_mean = np.log10(mean)
        ax.plot(f, log_mean, color=color, linewidth=2)
        ax.fill_between(f, np.log10(np.maximum(mean - ci, 1e-6)),
                        np.log10(mean + ci), color=color, alpha=0.25)
        handles.append(Line2D([0], [0], color=color, lw=2,
                              label=f"{pathology} (n={len(psds)})"))
    ax.set_xlim(0, 150)
    ax.set_ylim(-3, 4)
    ax.set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
    ax.set_title("PSD of cNS vs cNS+GBM", fontsize=10, fontweight='bold')
    ax.set_xlabel('Frequency (Hz)', fontsize=9)
    ax.set_ylabel('PSD\n(µV²/Hz)', fontsize=9)
    ax.legend(handles=handles, fontsize=8, loc='upper right', framealpha=0.7)
    ax.tick_params(labelsize=8)
    fig.tight_layout()

    pdf_path = os.path.join(OUTPUT_DIR, "PSD_GBM_vs_NeuronOnly.pdf")
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    print("=" * 60)
    print("ORGANOID PSD: GBM vs NEURON ONLY")
    print("=" * 60)

    recordings = discover_recordings()
    print(f"Found {len(recordings)} recordings\n")

    print("Loading / computing PSDs...")
    psd_data = load_all_psds(recordings)
    print(f"\nPSDs available for {len(psd_data)}/{len(recordings)} recordings\n")

    plot_gbm_vs_neuron(recordings, psd_data)

    print("\nDone!")
