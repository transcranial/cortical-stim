function [bestChan_bandPower, bestChan_tempCoher, bestChan_combined] = findBestChannel(signal_array, f_s, f_L, f_H, Nelectrodes)

%%%% Find channel with highest band power, highest temporal coherence, and
%%%% highest combined metric %%%%
%
% INPUTS
%   signal_array - input data as an array of signals
%   f_s - sampling frequency
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   Nelectrodes - number of channels in signal_array
%
% OUTPUTS
%   bestChan_bandPower - channel with highest power within frequency band
%   bestChan_tempCoher - channel with highest temporal coherence of
%       bandpass filtered signals
%   bestChan_combined - channel with highest combined metric (band power
%       and temporal coherence)
%%%%

freqs = f_L:0.01:f_H;
filterOrder_eval = 2;
ARorder_eval = 50;
for n_electrode=1:Nelectrodes
    signal = signal_array(n_electrode);
    Hs = spectrum.yulear(ARorder_eval);
    bandPower(n_electrode) = avgpower(psd(Hs,signal,'SpectrumType','onesided','Fs',f_s,'FreqPoints','User Defined','FrequencyVector',freqs));
    [b,a] = butter(filterOrder_eval,[f_L f_H]./(f_s/2),'bandpass');
    signal_filtered=filtfilt(b,a,signal);
    tempCoher(n_electrode) = temporalCoherence(signal_filtered,f_s);
end

[~, bestChan_bandPower] = max(bandPower);
[~, bestChan_tempCoher] = max(tempCoher);
[~, bandPower_sortedChan] = sort(bandPower,'descend');
[~, tempCoher_sortedChan] = sort(tempCoher,'descend');
for i=1:Nelectrodes
    x(i,1)=find(bandPower_sortedChan==i);
    y(i,1)=find(tempCoher_sortedChan==i);
end
[~, bestChan_combined] = min(x+y);

end