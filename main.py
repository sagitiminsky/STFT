# -*- coding: utf-8 -*-
"""HW2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1npPlJQmhi587E-d3Tbx2xCkZozZl_2Oe

# Imports
"""

import numpy as np
from scipy.io import wavfile, loadmat
from IPython.display import Audio
import matplotlib.pyplot as plt
from numpy.fft import fftshift
from STFT import STFT
from decode import *

DEBUG=1

"""# section 6"""
def section6():
    sample_rate, data = wavfile.read(r"./touchtone_1.wav")
    data = data / np.max(np.abs(data))
    Audio(data, rate=sample_rate)


"""# section 7"""
def section7():
    sample_rate, data = wavfile.read(r"./touchtone_1.wav")
    data = data / np.max(np.abs(data))

    plt.figure(figsize=(20, 3))
    freqs = np.fft.fftshift(np.fft.fftfreq(len(data), 1 / sample_rate))
    plt.plot(freqs, np.abs(np.fft.fftshift(np.fft.fft(data))))
    plt.xticks(np.arange(-2000, 2000, 500), np.arange(-2000, 2000, 500))
    plt.title("Frequency Domain", fontsize=16)
    plt.xlabel("Frequency [Hz]", fontsize=16)
    plt.ylabel("$|X^f(f)|$", fontsize=16)
    plt.grid()


"""# Section 8"""
def section8():
    sample_rate, data = wavfile.read(r"./touchtone_1.wav")
    data = data / np.max(np.abs(data))

    fig, ax = plt.subplots(2, 2, figsize=(16, 10))

    for row, col, win, hop, F in zip([0, 0, 1, 1], [0, 1, 0, 1],
                                     [512, 256, 1024, 2048], [512, 256, 1024, 2048],
                                     [64, 256, 256, 256]):
        X_stft = STFT(data, win, hop, F, sample_rate)
        tau = np.arange(X_stft.shape[1]) * hop / sample_rate
        freqs = np.fft.fftshift(np.fft.fftfreq(F, 1 / sample_rate))
        im = ax[row, col].pcolormesh(tau, freqs, fftshift(np.abs(X_stft), axes=0))
        ax[row, col].set_ylabel('f [Hz]', fontsize=16)
        ax[row, col].set_xlabel('$\\tau$ [sec]', fontsize=16)
        ax[row, col].set_title('win: ' + str(win) + '   hopSize: ' + str(hop) + '   F: ' + str(F), fontsize=16)
        fig.colorbar(im, ax=ax[row, col])

    fig.suptitle('| STFT(f, $\\tau$) |', fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])


"""# Section 9"""
def freqs2digit(freq1, freq2):
    digitsTable = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    highFreqs = np.array([1209, 1336, 1477])
    lowFreqs = np.array([697, 770, 852])

    minFreq = min(freq1, freq2)
    maxFreq = max(freq1, freq2)

    digit = digitsTable[np.abs(lowFreqs - minFreq).argmin(), np.abs(highFreqs - maxFreq).argmin()]

    return digit


def getCurrFreqs(stft_column, freqs, freq_dist_th, detected_th_param):
    allFreqs = np.array([697, 770, 852, 1209, 1336, 1477])
    stft_column = 20 * np.log10(stft_column + 1e-17)  # convert to DB
    detected_th = detected_th_param + np.median(stft_column)

    foundFreqsVal = []

    for search_freq in allFreqs:
        relevantFreqInds = (freqs > 0) & (np.abs(freqs - search_freq) < freq_dist_th)
        foundFreqsVal.append(stft_column[relevantFreqInds].mean())

    foundFreqsVal = np.array(foundFreqsVal)
    twoStrongetFreqsInds = foundFreqsVal.argsort()[-3:][::-1]

    freqs = []
    if foundFreqsVal[twoStrongetFreqsInds[1]] > detected_th:
        freqs = allFreqs[twoStrongetFreqsInds]

    return freqs





"""# Section 10"""
def section10():
    win = 256
    hop = 256
    F = 256
    freq_dist_th = 15
    detected_th_param = 8

    sample_rate, data = wavfile.read("./touchtone_1.wav")
    data = data / np.max(np.abs(data))
    t1_digits = decode(data, sample_rate, win, hop, F, freq_dist_th, detected_th_param)
    print('touchtone_1 digits: ' + str(t1_digits))


"""# Section 11"""
def section11():
    sample_rate, data = wavfile.read("./touchtone_2.wav")
    data = data / np.max(np.abs(data))
    digits_ground_truth = loadmat('./touchtone_2_sequence.mat')['sequence'].reshape(-1)

    fig, ax = plt.subplots(2, 3, figsize=(16, 10))
    win = 256
    hop = 256
    F = 256
    freq_dist_th = 15
    detected_th_param = 8

    err_rate = []
    sigma_arr = [0.05, 0.1, 0.25, 0.5, 1, 2] # [0.05, 0.1, 0.25, 0.5, 1, 2]

    # for detected_th_param in [5,6,7,8,9,10,11,12,13,14]:
    #     print("### detected_th_param:{} ###".format(detected_th_param))
    for ax, sigma in zip(fig.get_axes(), sigma_arr):
        noised_data = np.random.normal(0, sigma, len(data)) + data
        digits = decode(noised_data, sample_rate, win, hop, F, freq_dist_th, detected_th_param)
        cur_error_rate = 1 - np.sum([x == y for x, y in zip(np.array(digits), digits_ground_truth)]) / len(
            digits_ground_truth)

        if DEBUG==1:
            print('sigma:{} , digits len:{}, GT len:{}, current_error_rate:{} '.format(sigma, str(len(digits)),
                                                                                   str(len(digits_ground_truth)),
                                                                                   str(cur_error_rate)))
        err_rate.append(cur_error_rate)
        X_stft = STFT(noised_data, win, hop, F, sample_rate)
        tau = np.arange(X_stft.shape[1]) * hop / sample_rate
        freqs = np.fft.fftshift(np.fft.fftfreq(F, 1 / sample_rate))
        im = ax.pcolormesh(tau, freqs, fftshift(np.abs(X_stft), axes=0))
        ax.set_ylabel('f [Hz]', fontsize=16)
        ax.set_xlabel('$\\tau$ [sec]', fontsize=16)
        ax.set_title('sigma:{}'.format(sigma))
        fig.colorbar(im, ax=ax)

    plt.figure()
    plt.plot(sigma_arr, err_rate, '-*')
    plt.title('Error Rate')
    plt.xlabel('Sigma')
    plt.ylabel('Error rate')


"""# Section 12"""
def section12():
    win = 256
    hop = 256
    F = 256
    freq_dist_th = 15
    detected_th_param = 7

    sample_rate, data = wavfile.read("./touchtone_3.wav")
    data = data / np.max(np.abs(data))
    t3_digits = decode(data, sample_rate, win, hop, F, freq_dist_th, detected_th_param)
    print('touchtone_3 digits: ' + str(t3_digits))


if "__main__"==__name__:
    print("### section 6 ###")
    section6()
    print("### section 7 ###")
    section7()
    print("### section 8 ###")
    section8()
    print("### section 10 ###")
    section10()
    print("### section 11 ###")
    section11()
    print("### section 12 ###")
    section12()

    plt.show()