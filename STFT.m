function [S] = STFT(signal, win, hopSize, F,Fs )
%{
Input:
%%%%%%
    signal: input audio signal.
    win: the window to be used for STFT, if win is a scalar a Hamming widow of length win should be used.
    hopSize: time-hop between consecutive STFT time frames.
    F: if F is a scalar it should be the number of frequencies to be used by the FFT , in this case S should have F rows.
        If F is a vector it should contain the frequencies for which the STFT should be calculated at each time-step,
        in this case S should have the same number of rows as the number of elements in F.
    Fs: the signals sampling frequency.
Output:
%%%%%%%
    S: the STFT of the input signal. S should always have a number of columns equal to the number of time-steps in the signal.
%}

% check structure of F
if isscalar(F)
    S_rows = F;
    freq_vec = (0:(F-1))*2*pi/F;
    if F > length(signal) % pad signal if it's shorter than num of frequencies
        signal = [signal zeros(1, F-length(signal))];
    end
else % F is a vector
    if iscolumn(F) % convert to row if needed
        F = F';
    end
    freq_vec = F;
    if length(F) > length(signal) % pad signal if it's shorter than num of frequencies
        signal = [signal zeros(1, length(F)-length(signal))];
    end
    S_rows = length(F);
end

if isscalar(win)
    win_size = win;
    win = hamming(win_size); % returns a column
else
    if isrow(win) % convert to column if needed
        win = win';
    end
    win_size= length(win);
end

S_cols = floor((length(signal) - (win_size-hopSize))/hopSize);
S = zeros(S_rows, S_cols);

for ind_col = 1:S_cols

    % chop signal in relevant area and apply window on chopped signal:
    x = signal(((ind_col-1)*hopSize)+1 :((ind_col-1)*hopSize)+win_size);
    x = x .* win; 
    
    % perform DFT for required frequencies:
    phase_vec = 0:(length(x)-1);
    tmp_mat = freq_vec'*(-1i.*phase_vec);
    exp_mat = exp(tmp_mat);
    X_w = exp_mat*x;
    X_w = fftshift(X_w);
    
    S(:,ind_col) = X_w;
end

end