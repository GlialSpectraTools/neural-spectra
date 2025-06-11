import os
import numpy as np
import pandas as pd
from fooof import FOOOF
from scipy.interpolate import interp1d


def process_fooof(directory, output_dir, ID):
    
    results = []

    
    freq_file = os.path.join(directory, 'frequencies_psd.npy')
    if not os.path.exists(freq_file):
        raise FileNotFoundError("Shared frequency file 'frequencies_psd.npy' not found.")
    freqs = np.load(freq_file).flatten()
    if np.any(np.isnan(freqs)) or np.any(np.isinf(freqs)) or freqs.size == 0:
        raise ValueError("Frequency file contains invalid or empty values.")

    
    for filename in os.listdir(directory):
        if filename.endswith('_psd.npy'):
            
            parts = filename.split('_')
            electrode_name = parts[0]  

            
            if 'frequencies' in filename:
                continue

            
            modifying_number = parts[-2] if len(parts) > 1 else None  

            
            if modifying_number is None:
                raise ValueError(f"Cannot extract modifying number from filename {filename}")
            try:
                modifying_number = int(modifying_number)  
            except ValueError:
                raise ValueError(f"Modifying number '{modifying_number}' could not be converted to an integer.")

           
            psd = np.load(os.path.join(directory, filename)).flatten()
            if psd.size == 0 or np.any(np.isnan(psd)) or np.any(np.isinf(psd)):
                
                results.append({'Electrode': electrode_name, 'Slope': "empty", 'Offset': "empty", 'Modifying Number': modifying_number})
                print(f"Skipping {filename} due to empty or invalid values in PSD")
                continue

            
            mask = freqs != 120  
            freqs_filtered = freqs[mask]
            psd_filtered = psd[mask]

            interpolation = interp1d(freqs_filtered, psd_filtered, kind='linear', fill_value="extrapolate")
            freqs_interpolated = np.linspace(freqs_filtered.min(), freqs_filtered.max(), len(freqs_filtered))
            psd_interpolated = interpolation(freqs_interpolated)

            fm = FOOOF(peak_width_limits=[0.5, 12.0],
                       min_peak_height=0.05,
                       peak_threshold=2.0,
                       aperiodic_mode='fixed',
                       verbose=False)
            fm.fit(freqs_interpolated, psd_interpolated, [70, 150])

            slope = fm.get_params('aperiodic_params', 'exponent')
            offset = fm.get_params('aperiodic_params', 'offset')

            results.append({'Electrode': electrode_name, 'Slope': slope, 'Offset': offset, 'Modifying Number': modifying_number})

    df = pd.DataFrame(results)
    print("DataFrame before grouping:", df)  

    avg_params = df[df['Slope'] != "empty"].groupby('Modifying Number')[['Slope', 'Offset']].mean().reset_index()
    avg_params.columns = ['Modifying Number', 'Average Slope', 'Average Offset']

    df = pd.merge(df, avg_params, on='Modifying Number', how='left')

    output_file = os.path.join(output_dir, f'{ID}_fooof_results.csv')

    df.to_csv(output_file, index=False)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    directory = #'directory to psd folder from script #1'
    output_dir = #'output directory of choice'
    ID = #'PtID_correctVincorrect_psd_output'
    process_fooof(directory, output_dir, ID)
