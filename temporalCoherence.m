function [tau_coherence] = temporalCoherence(signal, f_s)

%%%% Temporal Coherence of signal %%%%

c=xcorr(signal);
c_sym = c(length(signal):(2*length(signal)-1));
c_max = c(length(signal));
index=find(c_sym>0.5*c_max,1,'last');

tau_coherence = index/f_s;

end