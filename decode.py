import numpy as np
from numpy.fft import fftshift
from STFT import STFT
from main import getCurrFreqs
from main import freqs2digit

def decode(signal, sample_rate=4096, win=256, hop=256, F=256, freq_dist_th=30, detected_th_param=8):
    freqs = np.fft.fftshift(np.fft.fftfreq(F, 1 / sample_rate))

    digits = []
    waitForQuite = False
    for col_idx in range(len(signal) // win):
        curr_col_stft = STFT(signal[(col_idx * hop):(col_idx * hop + win)], win, hop, F, sample_rate)
        curr_col_stft = fftshift(np.abs(curr_col_stft), axes=0)
        chosen_freqs = getCurrFreqs(curr_col_stft, freqs, freq_dist_th, detected_th_param)
        if len(chosen_freqs) > 0 and not waitForQuite:
            digits.append(freqs2digit(chosen_freqs[0], chosen_freqs[1]))
            waitForQuite = True
        elif len(chosen_freqs) == 0:
            waitForQuite = False

    return digits