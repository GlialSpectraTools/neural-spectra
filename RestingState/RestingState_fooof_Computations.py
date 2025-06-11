import os
import numpy as np
import pandas as pd
from fooof import FOOOF
from scipy.interpolate import interp1d


def process_fooof(directory, output_dir, ID):

    results = []

 
    for filename in os.listdir(directory):
        if filename.endswith('_frequencies.npy'):
           
            electrode_name = filename.split('_')[0]  

           
            freqs = np.load(os.path.join(directory, filename))
            psd_filename = filename.replace('_frequencies.npy', '_psd.npy')
            psd = np.load(os.path.join(directory, psd_filename))

         
            freqs = freqs.flatten()
            psd = psd.flatten()

        
            if freqs.ndim != 1 or psd.ndim != 1:
                raise ValueError(f"freqs and psd for electrode {electrode_name} should both be 1D arrays")
            if np.any(np.isnan(freqs)) or np.any(np.isnan(psd)):
                raise ValueError(f"NaN values found in electrode {electrode_name}")
            if np.any(np.isinf(freqs)) or np.any(np.isinf(psd)):
                raise ValueError(f"Infinite values found in electrode {electrode_name}")
            if freqs.size == 0 or psd.size == 0:
                raise ValueError(f"Empty array for electrode {electrode_name}")

           
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
            fm.fit(freqs_interpolated, psd_interpolated, [1, 150])  # Change range to frequency of interest

          
            slope = fm.get_params('aperiodic_params', 'exponent')
            intercept = fm.get_params('aperiodic_params', 'offset')


            results.append({'Electrode': electrode_name, 'Slope': slope, 'Intercept': intercept})

  
    df = pd.DataFrame(results)

   
    output_file = os.path.join(output_dir, f'{ID}_fooof_results.csv')

    
    df.to_csv(output_file, index=False)

    print(f"Results saved to {output_file}")


# Usage example
if __name__ == "__main__":
    directory = #'Paste directory of PSD folder from Welch's method script'  
    output_dir = #'Paste output directory'
    ID = '#'  # Specify your ID here
    process_fooof(directory, output_dir, ID)
