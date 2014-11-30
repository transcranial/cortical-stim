function [signal_segment_filtered] = zeroPhaseFilter(signal_segment, order, type, f_L, f_H, f_s)

%%%% Zero-Phase Bandpass Filter %%%%
%
% INPUTS
%   signal_segment - input signal segment
%   order - order of filter
%   type - filter type (butterworth, chebyshev1, chebyshev2, elliptic,
%       bessel)
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   f_s - sampling frequency
%
% OUTPUTS
%   signal_segment_filtered - output filtered signal segment
%%%%

switch type
    case 'butterworth'
        [b,a] = butter(order,[f_L f_H]./(f_s/2),'bandpass');
    case 'chebyshev1'
        [b,a] = cheby1(order,0.1,[f_L f_H]./(f_s/2),'bandpass');
    case 'chebyshev2'
        [b,a] = cheby2(order,60,[f_L f_H]./(f_s/2),'bandpass');
    case 'elliptic'
        [b,a] = ellip(order,1,60,[f_L f_H]./(f_s/2),'bandpass');
    case 'bessel'
        [b,a] = besself(order,[f_L f_H]./(f_s/2),'bandpass');
end

signal_segment_filtered = filtfilt(b,a,signal_segment);

end