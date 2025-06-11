import mne
import numpy as np
import os
from scipy.signal import welch, periodogram
from scipy import fft as sp_fft
from scipy import signal
import pandas as pd

def simplepwr(s,fs,npts = None): 
    N = len(s)
    if npts is None:
        npts = len(s) 

    nyquist = fs/2
    freqs = np.linspace(0,int(nyquist),int((npts/2) + 1))

    
    s = s - np.mean(s)

    
    windowed_data = signal.windows.hann(N)*s
    raw_fft = sp_fft.rfft(windowed_data, npts)
    pwr_init = np.conjugate(raw_fft)*raw_fft
    scale = 1/((N/2)**2)
    pwr_scaled = pwr_init*scale

    if npts % 2:
        pwr_scaled[..., 1:] *= 2
    else:
        
        pwr_scaled[..., 1:-1] *= 2

    return freqs, pwr_scaled.real


def load_eeg_set(file_path):
    
    return mne.io.read_raw_eeglab(file_path, preload=True)


def notch_filter(eeg_data, freqs=[60, 120, 180, 240]):
    
    eeg_data.notch_filter(freqs=freqs, fir_design='firwin')
    return eeg_data


def calculate_welch_psd(eeg_data, sfreq, nperseg=None):
    
    data = eeg_data.get_data()
    n_channels = data.shape[0]

    
    nperseg = nperseg or int(sfreq * 0.5)
    data_length = data.shape[1]

    
    nperseg = min(nperseg, data_length)

    
    noverlap = max(nperseg // 2, 1)

    psd_list = []
    frequencies = None
   
    for i in range(n_channels):
        channel_data = data[i, :]
        
        freqs,psd = periodogram(channel_data, fs=sfreq, window='hann', scaling='spectrum')

        if frequencies is None:
            frequencies = freqs

        psd_list.append(psd)

    return frequencies, psd_list


def calculate_power_fft(eeg_data, sfreq):
    
    data = eeg_data.get_data()
    n_channels = data.shape[0]

    
    spectrum_list = []
    frequencies = None
    
    for i in range(n_channels):
        channel_data = data[i, :]
        freqs, spectrum = simplepwr(channel_data, sfreq) 
        if frequencies is None:
            frequencies = freqs
        spectrum_list.append(spectrum)

    return frequencies, spectrum_list


def save_frequencies_and_psd(output_dir, channel_names, frequencies, psd_averages, trigger_value, save_freqs = 1):
    
    os.makedirs(output_dir, exist_ok=True)

    if save_freqs == 1:
       
        np.save(os.path.join(output_dir, 'frequencies_psd.npy'), frequencies)

    for i, channel_name in enumerate(channel_names):
        
        np.save(os.path.join(output_dir, f'{channel_name}_trigger_{trigger_value}_psd.npy'), psd_averages[i])


def process_ecog_set(file_path, output_folder_name, trigger_values, post_trigger_duration=0.5):
    
    raw = load_eeg_set(file_path)
    raw = notch_filter(raw)

    
    events, _ = mne.events_from_annotations(raw)
    events = events[1:]

    
    if len(trigger_values) != len(events):
        raise ValueError("Number of trigger values must match the number of events (excluding the first).")

    
    sfreq = raw.info['sfreq']
    post_trigger_samples = int(post_trigger_duration * sfreq)

    
    output_dir = os.path.join(os.path.dirname(file_path), output_folder_name)
    os.makedirs(output_dir, exist_ok=True)

    
    psd_zeros = [[] for _ in range(len(raw.ch_names))]
    psd_ones = [[] for _ in range(len(raw.ch_names))]
    frequencies = None  

    
    for idx, event in enumerate(events):
        trigger_sample = event[0]
        trigger_value = trigger_values[idx]

        
        if trigger_sample + post_trigger_samples <= len(raw.times):
            segment_data = raw[:, trigger_sample:trigger_sample + post_trigger_samples][0]

            
            temp_info = mne.create_info(ch_names=raw.ch_names, sfreq=sfreq, ch_types='eeg')
            temp_raw = mne.io.RawArray(segment_data, temp_info)

            
            freqs, psd = calculate_welch_psd(temp_raw, sfreq, nperseg=int(sfreq * 0.5))
            

            
            if frequencies is None:
                frequencies = freqs

            
            if trigger_value == 0:
                for i, p in enumerate(psd):
                    psd_zeros[i].append(p)
            elif trigger_value == 1:
                for i, p in enumerate(psd):
                    psd_ones[i].append(p)

            print(f"Processed event {idx + 1}/{len(events)} with trigger {trigger_value}")

    
    averaged_psd_zeros = [np.mean(p, axis=0) for p in psd_zeros]
    averaged_psd_ones = [np.mean(p, axis=0) for p in psd_ones]

    
    save_frequencies_and_psd(output_dir, raw.ch_names, frequencies, averaged_psd_zeros, trigger_value=0, save_freqs = 1)
    save_frequencies_and_psd(output_dir, raw.ch_names, frequencies, averaged_psd_ones, trigger_value=1, save_freqs = 0)


#file_path = '.set input directory'
output_folder_name = #'patientID_correctVincorrect_psd_output' 
#trigger values: 0 = correct, 1 = incorrect
trigger_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] #Change to x amount of stimuli for PN and x amount of stimuli for AN


process_ecog_set(file_path, output_folder_name, trigger_values)
