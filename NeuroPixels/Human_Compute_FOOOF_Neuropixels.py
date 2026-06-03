#!/usr/bin/env python3
"""
FOOOF (1/f) Computation for Neuropixels LFP Data - Simplified Version
======================================================================
Computes FOOOF slope (exponent) for each channel from Neuropixels LFP data.
Uses direct PSD computation over the entire interest window for speed.
"""

import os
import numpy as np
import pandas as pd
from fooof import FOOOF
from scipy.interpolate import interp1d
from scipy.signal import welch
from tqdm import tqdm

# ==============================================================
# Configuration — set BASE_DIR to the root of the project
# ==============================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

metadata_path = os.path.join(BASE_DIR, "metadata", "np_meta_sheet.xlsx")
data_base_path = os.path.join(BASE_DIR, "raw_lfp")
output_dir = os.path.join(BASE_DIR, "patient_csvs")
os.makedirs(output_dir, exist_ok=True)

fs = 400  # Sampling rate in Hz

# ==============================================================
# Helper Functions
# ==============================================================

def remove_line_noise(freqs, psd):
    """Removes 60Hz line noise and harmonics from PSD."""
    clean_psd = psd.copy()
    for notch_freq in [60, 120, 180]:
        notch_idx = np.where((freqs >= notch_freq - 3) & (freqs <= notch_freq + 3))[0]
        if len(notch_idx) > 0 and notch_idx[0] > 0 and notch_idx[-1] < len(psd) - 1:
            left_idx = max(0, notch_idx[0] - 5)
            right_idx = min(len(freqs) - 1, notch_idx[-1] + 5)
            clean_psd[notch_idx] = np.interp(
                notch_idx, [left_idx, right_idx], [clean_psd[left_idx], clean_psd[right_idx]]
            )
    return clean_psd

def process_fooof_for_patient(signal_file, distance_file, output_dir, patient_id, block,
                               interest_start_sec, interest_end_sec, fs=400):
    """
    Process FOOOF for a single patient/block/probe.
    Computes PSD over entire interest window, then fits FOOOF per channel.
    """
    print(f"\nProcessing {patient_id}_{block}...")
    
    # Load LFP data
    lfp = np.load(signal_file, mmap_mode='r')
    if lfp.ndim == 2 and lfp.shape[0] < lfp.shape[1]:
        lfp = lfp.T  # Ensure (time, channels)
    
    n_samples, n_channels = lfp.shape
    total_duration = n_samples / fs
    
    # Handle 'inf' end time
    if interest_end_sec == 'inf' or interest_end_sec > total_duration:
        interest_end_sec = total_duration
    
    print(f"  LFP: {n_channels} channels, {total_duration:.1f}s total")
    print(f"  Interest window: {interest_start_sec}s to {interest_end_sec:.1f}s")
    
    # Extract interest window
    start_idx = int(interest_start_sec * fs)
    end_idx = int(interest_end_sec * fs)
    lfp_interest = np.array(lfp[start_idx:end_idx, :])  # Copy to avoid mmap issues
    
    print(f"  Extracted {lfp_interest.shape[0] / fs:.1f}s of data")
    
    # Load distances
    distances = np.load(distance_file) if os.path.exists(distance_file) else None
    
    # Fix for patients with channel indices instead of µm distances
    # If max distance < 1000, assume it's channel indices and convert to µm (20µm spacing)
    if distances is not None and np.max(distances) < 1000:
        print(f"  Warning: Converting channel indices to µm (20µm spacing)")
        distances = distances * 20
    
    # Compute PSD for all channels at once using Welch's method
    # Using 3-second windows with 50% overlap to match ECoG paper methodology
    print("  Computing PSDs...")
    nperseg = min(fs * 3, lfp_interest.shape[0])  # 3-second windows (matching ECoG paper)
    f, Pxx = welch(lfp_interest, fs=fs, nperseg=nperseg, noverlap=nperseg//2, axis=0)
    # Pxx shape: (frequencies, channels)
    
    # Fit FOOOF for each channel
    print(f"  Fitting FOOOF on {n_channels} channels...")
    results = []
    
    for ch_idx in tqdm(range(n_channels), desc=f"  {patient_id}_{block}", leave=False):
        psd = Pxx[:, ch_idx]
        
        # Skip bad channels
        if np.any(np.isnan(psd)) or np.all(psd <= 0):
            continue
        
        # Remove line noise
        psd_clean = remove_line_noise(f, psd)
        psd_clean = np.maximum(psd_clean, 1e-15)  # Ensure positive
        
        try:
            fm = FOOOF(peak_width_limits=[0.5, 12.0],
                      min_peak_height=0.05,
                      peak_threshold=2.0,
                      aperiodic_mode='fixed',
                      max_n_peaks=6,
                      verbose=False)
            fm.fit(f, psd_clean, [70, 150])  # Fit 70-150 Hz range (high gamma)
            
            slope = fm.get_params('aperiodic_params', 'exponent')
            distance_to_tip = distances[ch_idx] if distances is not None else np.nan
            
            results.append({
                'Channel': ch_idx,
                'Distance_to_Tip': distance_to_tip,
                'Slope': slope
            })
        except:
            continue
    
    if not results:
        print(f"  Warning: No valid results for {patient_id}_{block}")
        return None
    
    # Save CSV
    df = pd.DataFrame(results)
    imec_id = os.path.basename(signal_file).split('_imec')[1].split('_')[0] if '_imec' in signal_file else '0'
    output_file = os.path.join(output_dir, f'{patient_id}_{block}_imec{imec_id}_fooof_slope_by_channel.csv')
    df.to_csv(output_file, index=False)
    
    print(f"  Saved {len(results)} channels to {os.path.basename(output_file)}")
    return df

# ==============================================================
# MAIN
# ==============================================================

if __name__ == "__main__":
    print("="*60)
    print("FOOOF COMPUTATION FOR NEUROPIXELS")
    print("="*60)
    
    # Load metadata
    df_metadata = pd.read_excel(metadata_path, sheet_name='Recordings')
    available_patients = [d for d in os.listdir(data_base_path) 
                         if d.startswith('NP') and os.path.isdir(os.path.join(data_base_path, d))]
    print(f"Found {len(available_patients)} patients")
    
    # Parse interest windows from metadata
    interest_windows = {}
    for _, row in df_metadata.iterrows():
        subject, block, sort_time = row['Subject'], row['Block'], row['Sort time']
        if pd.notna(sort_time) and subject in available_patients:
            try:
                if isinstance(sort_time, str):
                    clean = sort_time.replace('[', '').replace(']', '').replace(' ', '')
                    if ',' in clean:
                        start, end = clean.split(',')
                        end = 'inf' if end.lower() == 'inf' else int(end)
                        interest_windows[f"{subject}_{block}"] = [int(start), end]
            except:
                pass
    
    print(f"Found {len(interest_windows)} interest windows")
    
    # Process each patient
    processed = 0
    for subject in sorted(available_patients):
        subject_dir = os.path.join(data_base_path, subject)
        
        # Find matching blocks
        for block_key, (start_sec, end_sec) in interest_windows.items():
            if not block_key.startswith(subject + '_'):
                continue
            
            block = block_key.split('_')[1]
            
            # Find LFP files for this block
            lfp_files = [f for f in os.listdir(subject_dir) 
                        if f.endswith('_lfp_filtered_400Hz.npy') and f.startswith(f"{subject}_{block}")]
            
            for lfp_file in lfp_files:
                imec_id = lfp_file.split('_imec')[1].split('_')[0] if '_imec' in lfp_file else '0'
                dist_file = f"{subject}_{block}_g0_imec{imec_id}_channel_distances_to_tip.npy"
                
                signal_path = os.path.join(subject_dir, lfp_file)
                dist_path = os.path.join(subject_dir, dist_file)
                
                if not os.path.exists(dist_path):
                    print(f"Skipping {lfp_file} - no distance file")
                    continue
                
                try:
                    result = process_fooof_for_patient(
                        signal_path, dist_path, output_dir,
                        subject, block, start_sec, end_sec, fs
                    )
                    if result is not None:
                        processed += 1
                except Exception as e:
                    print(f"Error processing {block_key}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Done! Processed {processed} files")
    print(f"Results saved to: {output_dir}")
    print("="*60)
