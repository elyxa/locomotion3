function filtered_channels = temporal_filter(channels, Fs)
% temporal filtering of the signal 
% (girls' code replaced by Andr??s's butter filter code...)

D = designfilt('lowpassfir','FilterOrder',4,'PassbandFrequency',4.5, 'StopbandFrequency', 5, 'SampleRate',Fs);

filtered_channels = filter(D,channels);
filtered_channels = filtered_channels;

%Cite frequencies: CAN SACRAL MARKER APPROXIMATE CENTER OF MASS DURING GAIT AND SLIP-FALL RECOVERY AMONG COMMUNITY-DWELLING OLDER ADULTS?

end
