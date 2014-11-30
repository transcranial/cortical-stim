function [stimDelivered, t_delay] = stimRT(signal_segment, f_s, f_L, f_H, fstepsize, ARorder, lambda, filterOrder, filterType, ARmethod, t_start, t_stop, targetPhase)

%%%% Phase-Locked Stimulation (Real Time) %%%%
%
% INPUTS
%   signal_segment - input signal segment (single channel input)
%   f_s - sampling frequency
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   fstepsize - frequency step size for passband optimization
%   ARorder - order of autoregressive model
%   lambda - threshold fraction of power in optimized frequency band
%   filterOrder - order of bandpass filter
%   filterType - type of bandpass filter
%   ARmethod - method used to build AR model
%   t_start - start time of portion of signal used for AR model, in sec
%   t_stop - stop time of portion of signal used for AR model, in sec
%   targetPhase - (0 corresponds to peak, 180 corresponds to trough)
%
% OUTPUTS
%   stimDelivered - boolean for whether or not stimulation is delivered
%   t_delay - time delay for stimulation output, in sec (only valid if
%       stimDelivered = 1)
%%%%

%%%% Example Parameter Initialization %%%%
% f_s = 500;
% f_L = 4;
% f_H = 9;
% fstepsize = 0.05;
% ARorder = optParams(i,1);
% lambda = optParams(i,2);
% filterOrder = optParams(i,3);
% filterType = filters{optParams(i,4)};
% ARmethod = 'yw';
% t_stop = optParams(i,5);
% t_start = segmentLength - t_stop;
% targetPhase = 0;
% filters = {'butterworth', 'chebyshev1', 'chebyshev2', 'elliptic', 'bessel'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frequencies = f_L:fstepsize:f_H;
Pxx = pyulear(signal_segment, ARorder, frequencies, f_s); 

stimDelivered = 0;
t_delay = 0;
if ~isempty(findpeaks(Pxx))
    [fhat_L, fhat_H] = freqBandOpt(signal_segment, f_L, f_H, f_s, fstepsize, ARorder, lambda);
    [signal_segment_filtered] = zeroPhaseFilter(signal_segment, filterOrder, filterType, fhat_L, fhat_H, f_s);
    [signal_predicted, index0] = linearPredAR(signal_segment_filtered, ARorder, ARmethod, f_s, t_start, t_stop);
    [phi_inst, f_inst] = instPhaseFreq(signal_predicted, index0, f_s);
    t_delay = (1/f_inst) * mod(targetPhase-phi_inst+360,360) / 360;

    if t_delay>0 && t_delay<1/f_L        
        stimDelivered = 1;
    else
        stimDelivered = 0;
    end
end

end

