#!/usr/bin/env python3
"""
FOOOF (1/f) Computation for Organoid Neuropixels LFP Data
==========================================================
Loads organoid NWB files via pynapple, restricts to the first 3 minutes,
downsamples from ~2500 Hz to 400 Hz, then computes FOOOF aperiodic exponent
per channel over the 70-150 Hz (high gamma) range.
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import pynapple as nap
from pynwb import NWBHDF5IO
from fooof import FOOOF
from scipy.signal import welch, resample_poly
from tqdm import tqdm

# ==============================================================
# Configuration
# ==============================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NWB_BASE_PATH = os.path.join(BASE_DIR, "raw_nwb")
OUTPUT_DIR = os.path.join(BASE_DIR, "organoid_csvs", "recording_csvs")
COMBINED_CSV = os.path.join(BASE_DIR, "organoid_csvs", "organoid_oneoverf_all_channels.csv")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_FS = 400
LFP_DURATION_SEC = 180  # first 3 minutes

SAMPLE_METADATA = {
    "HCN-1":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "HCN-2":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "HCN-3":    {"group": "neuron_only", "pathology": "Neuron Only"},
    "SF0683-1": {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF0683-4": {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF0683-6": {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF0683-7": {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF496-1":  {"group": "tumor", "pathology": "Grade IV Astro"},
    "SF496-3":  {"group": "tumor", "pathology": "Grade IV Astro"},
    "SF702-1":  {"group": "tumor", "pathology": "Grade II Astro"},
    "SF702-3":  {"group": "tumor", "pathology": "Grade II Astro"},
    "SF718-1":  {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF718-2":  {"group": "tumor", "pathology": "Grade II Oligo"},
    "SF722-1":  {"group": "tumor", "pathology": "Grade III Oligo"},
    "SF722-2":  {"group": "tumor", "pathology": "Grade III Oligo"},
    "SF722-3":  {"group": "tumor", "pathology": "Grade III Oligo"},
    "SF722-4":  {"group": "tumor", "pathology": "Grade III Oligo"},
    "SF622-1":  {"group": "tumor", "pathology": "Grade IV Astro"},
    "SF622-2":  {"group": "tumor", "pathology": "Grade IV Astro"},
    "SF725-1":  {"group": "tumor", "pathology": "GBM"},
    "SF725-2":  {"group": "tumor", "pathology": "GBM"},
    "SF725-3":  {"group": "tumor", "pathology": "GBM"},
    "SF725-4":  {"group": "tumor", "pathology": "GBM"},
}

# ==============================================================
# Helper Functions
# ==============================================================

def remove_line_noise(freqs, psd):
    """Removes 60 Hz line noise and harmonics from PSD via interpolation."""
    clean_psd = psd.copy()
    for notch_freq in [60, 120, 180]:
        notch_idx = np.where((freqs >= notch_freq - 3) & (freqs <= notch_freq + 3))[0]
        if len(notch_idx) > 0 and notch_idx[0] > 0 and notch_idx[-1] < len(psd) - 1:
            left_idx = max(0, notch_idx[0] - 5)
            right_idx = min(len(freqs) - 1, notch_idx[-1] + 5)
            clean_psd[notch_idx] = np.interp(
                notch_idx, [left_idx, right_idx],
                [clean_psd[left_idx], clean_psd[right_idx]]
            )
    return clean_psd


def discover_recordings():
    """Find all non-plx organoid recordings on disk that match SAMPLE_METADATA."""
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
            else:
                print(f"  Warning: NWB file missing in {dirname}")
        if not dirs:
            print(f"  Warning: No directories found for sample {sample}")
    return recordings


def process_recording(rec, target_fs=TARGET_FS, duration_sec=LFP_DURATION_SEC):
    """
    Load NWB, restrict to first `duration_sec` seconds, downsample to
    `target_fs` Hz, compute PSD + FOOOF per channel.
    """
    recording_name = rec["recording"]
    print(f"\nProcessing {recording_name}...")

    # Read ADC-to-volts conversion factor from NWB metadata
    io = NWBHDF5IO(rec["nwb_path"], 'r')
    nwbfile = io.read()
    conversion = nwbfile.acquisition['ElectricalSeriesLF'].conversion
    io.close()

    data = nap.load_file(rec["nwb_path"])
    lfp_raw = data["ElectricalSeriesLF"]
    native_fs = lfp_raw.rate
    total_dur = lfp_raw.time_support.end[0] - lfp_raw.time_support.start[0]
    print(f"  Native rate: {native_fs:.1f} Hz, {lfp_raw.shape[1]} channels, {total_dur:.1f}s total")
    print(f"  ADC conversion: {conversion:.6e} V/count")

    end_time = min(lfp_raw.time_support.start[0] + duration_sec,
                   lfp_raw.time_support.end[0])
    lfp_restricted = lfp_raw.restrict(
        nap.IntervalSet(start=lfp_raw.time_support.start[0], end=end_time)
    )
    # Convert from raw ADC int16 to microvolts
    lfp_np = np.asarray(lfp_restricted.values, dtype=np.float64) * conversion * 1e6
    actual_dur = lfp_np.shape[0] / native_fs
    print(f"  Restricted to {actual_dur:.1f}s ({lfp_np.shape[0]} samples), converted to µV")

    # Downsample to target_fs using rational resampling
    # native_fs ≈ 2500, target 400 -> up=4, down=25
    up, down = 4, 25
    effective_target = native_fs * up / down
    print(f"  Resampling {native_fs:.1f} -> {effective_target:.1f} Hz (up={up}, down={down})")
    lfp_ds = resample_poly(lfp_np, up, down, axis=0)
    n_samples, n_channels = lfp_ds.shape
    print(f"  Downsampled: {n_samples} samples, {n_channels} channels")

    # Welch PSD — 3-second windows, 50% overlap
    nperseg = min(target_fs * 3, n_samples)
    f, Pxx = welch(lfp_ds, fs=target_fs, nperseg=nperseg, noverlap=nperseg // 2, axis=0)

    # FOOOF per channel
    print(f"  Fitting FOOOF on {n_channels} channels...")
    results = []
    for ch_idx in tqdm(range(n_channels), desc=f"  {recording_name}", leave=False):
        psd = Pxx[:, ch_idx]
        if np.any(np.isnan(psd)) or np.all(psd <= 0):
            continue

        psd_clean = remove_line_noise(f, psd)
        psd_clean = np.maximum(psd_clean, 1e-15)

        try:
            fm = FOOOF(
                peak_width_limits=[0.5, 12.0],
                min_peak_height=0.05,
                peak_threshold=2.0,
                aperiodic_mode='fixed',
                max_n_peaks=6,
                verbose=False,
            )
            fm.fit(f, psd_clean, [70, 150])
            slope = fm.get_params('aperiodic_params', 'exponent')
            results.append({
                "Channel": ch_idx,
                "Slope": slope,
                "Recording": recording_name,
                "Sample": rec["sample"],
                "Group": rec["group"],
                "Pathology": rec["pathology"],
            })
        except Exception:
            continue

    if not results:
        print(f"  Warning: No valid FOOOF results for {recording_name}")
        return None

    df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, f"{recording_name}_fooof_slope.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved {len(results)} channels -> {os.path.basename(out_path)}")
    return df


# ==============================================================
# MAIN
# ==============================================================

def combine_csvs():
    """Combine all per-recording CSVs into a single summary CSV."""
    csv_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "*_fooof_slope.csv")))
    if not csv_files:
        print("No per-recording CSVs found to combine.")
        return
    dfs = [pd.read_csv(f) for f in csv_files]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(COMBINED_CSV, index=False)
    print(f"Combined {len(csv_files)} files ({len(combined)} rows) -> {COMBINED_CSV}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organoid 1/f FOOOF computation")
    parser.add_argument("--recording", type=str, default=None,
                        help="Process a single recording, e.g. NPorganoid_HCN-3_B1")
    parser.add_argument("--combine-only", action="store_true",
                        help="Only combine existing per-recording CSVs into summary CSV")
    args = parser.parse_args()

    print("=" * 60)
    print("ORGANOID 1/f (FOOOF) COMPUTATION")
    print("=" * 60)

    if args.combine_only:
        combine_csvs()
        sys.exit(0)

    recordings = discover_recordings()
    print(f"\nFound {len(recordings)} recordings across "
          f"{len(set(r['sample'] for r in recordings))} samples")

    if args.recording:
        matches = [r for r in recordings if r["recording"] == args.recording]
        if not matches:
            print(f"ERROR: Recording '{args.recording}' not found on disk.")
            sys.exit(1)
        recordings = matches

    all_dfs = []
    processed = 0
    for rec in recordings:
        try:
            df = process_recording(rec)
            if df is not None:
                all_dfs.append(df)
                processed += 1
        except Exception as e:
            print(f"  ERROR processing {rec['recording']}: {e}")

    if not args.recording:
        if all_dfs:
            combined = pd.concat(all_dfs, ignore_index=True)
            combined.to_csv(COMBINED_CSV, index=False)
            print(f"\nCombined CSV: {len(combined)} rows -> {COMBINED_CSV}")
        else:
            print("\nNo data to combine.")

    print(f"\n{'=' * 60}")
    print(f"Done! Processed {processed}/{len(recordings)} recordings")
    print("=" * 60)
