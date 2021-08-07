import numpy as np


"""# section 5"""
def STFT(signal, win, hopSize, F, Fs):
    if not hasattr(win, "__len__"):
        win = np.hamming(win)
    if not hasattr(F, "__len__"):
        F = 2*np.pi*np.arange(F)/F


    stft = []
    startIdx = 0
    while startIdx + len(win) <= len(signal):
        t = np.arange(startIdx, startIdx + len(win))
        e = np.exp(-1j * t.reshape(1, -1) * F.reshape(-1, 1))
        currDFT = np.sum(signal[startIdx:(startIdx + len(win))]*win*e, 1)
        stft.append(np.abs(currDFT))
        startIdx += hopSize

    stft = np.stack(stft).T
    return stft