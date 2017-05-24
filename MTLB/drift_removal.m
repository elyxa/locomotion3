function new_signal = drift_removal(signal)

len = length(signal);
detrended_signal = detrend(signal);
diff = signal(1) - detrended_signal(1);
new_signal = detrended_signal + diff*ones(1,len)

end