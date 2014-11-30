function [fhat_L, fhat_H] = freqBandOpt(signal_segment, f_L, f_H, f_s, fstepsize, ARorder, lambda)

%%%% Frequency Band Optimization %%%%
%
% INPUTS
%   signal_segment - input signal segment
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   f_s - sampling frequency
%   fstepsize - frequency step size for optimization
%   ARorder - order of autoregressive spectral estimation
%   lambda - threshold fraction of power in optimized frequency band
%
% OUTPUTS
%   fhat_L - optimized passband lower frequency limit
%   fhat_H - optimized passband higher frequency limit
%%%%

fhat_L = f_L;
fhat_H = f_H;

frequencies = f_L:fstepsize:f_H;

Pxx = pyulear(signal_segment, ARorder, frequencies, f_s);

[~,index_f_L] = min(abs(frequencies-f_L));
[~,index_f_H] = min(abs(frequencies-f_H));
[~,index_fhat_L] = min(abs(frequencies-fhat_L));
[~,index_fhat_H] = min(abs(frequencies-fhat_H));

while sum(Pxx(index_fhat_L:index_fhat_H)) > lambda*sum(Pxx(index_f_L:index_f_H))
    if Pxx(index_fhat_L) < Pxx(index_fhat_H)
        fhat_L = fhat_L + fstepsize;
        [~,index_fhat_L] = min(abs(frequencies-fhat_L));
    else
        fhat_H = fhat_H - fstepsize;
        [~,index_fhat_H] = min(abs(frequencies-fhat_H));
    end
end

end