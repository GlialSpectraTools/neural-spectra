import numpy as np
import mne 
from scipy.signal import welch, detrend
import os


# Load the .set EEG file using MNE
def load_eeg_set(file_path):
    return mne.io.read_raw_eeglab(file_path, preload=True)


def notch_filter(eeg_data, freqs=[60, 120, 180, 240]):
    eeg_data.notch_filter(freqs=freqs)
    return eeg_data


def calculate_welch_psd(eeg_data, sfreq, fmin=1, fmax=150, nperseg=None): #CHANGE FREQUENCIES HERE IF NEEDED
    data = eeg_data.get_data()
    n_channels = data.shape[0]

    nperseg = sfreq * 3  

    frequencies_list = []
    psd_list = []


    for i in range(n_channels):
        channel_data = data[i, :]
        freqs, psd = welch(channel_data, fs=sfreq, nperseg=nperseg, noverlap=nperseg // 2)

        frequencies_list.append(freqs)
        psd_list.append(psd)

    return frequencies_list, psd_list


def save_channel_data(output_dir, channel_names, frequencies, psd):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, channel_name in enumerate(channel_names):
        np.save(os.path.join(output_dir, f'{channel_name}_frequencies.npy'), frequencies[i])
        np.save(os.path.join(output_dir, f'{channel_name}_psd.npy'), psd[i])


def main(file_path, file_id):
    
    eeg_data = load_eeg_set(file_path)

 
    sfreq = eeg_data.info['sfreq']
    print(f"Sampling frequency: {sfreq} Hz")


    eeg_data = notch_filter(eeg_data) 

  
    channel_names = eeg_data.ch_names

  
    frequencies, psd = calculate_welch_psd(eeg_data, sfreq)

  
    output_dir = f'{file_id}_output_psd'

   
    save_channel_data(output_dir, channel_names, frequencies, psd)

    print(f"Frequencies and power spectra saved in '{output_dir}' for each channel.")


# Usage
if __name__ == "__main__":
    file_id = #'Insert file name'
    file_path = #f'insert .set file path'  # Replace with the path to your .set file
    main(file_path, file_id)
