function [phi_inst, f_inst] = instPhaseFreq(signal_segment, index0, f_s)

%%%% Instantaneous Phase and Frequency %%%%
%
% INPUTS
%   signal_segment - input signal segment
%   index0 - index of signal_segment where t=0
%   f_s - sampling frequency
%
% OUTPUTS
%   phi_inst - instantaneous phase, in degrees
%   f_inst - instantaneous frequency, in Hz
%%%%

H_segment = hilbert(signal_segment);
phi_segment = angle(H_segment);
f_segment = f_s.*diff(unwrap(phi_segment))./(2*pi);

phi_inst = phi_segment(index0).*180./pi;
f_inst = f_segment(index0);

end