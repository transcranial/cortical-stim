function [indexStimTimes, stimFreqs, stimPhases, phaseLockingMetrics] = stimSim(data, f_s, segmentLength, f_L, f_H, fstepsize, ARorder, lambda, filterOrder, filterType, ARorder_eval, filterOrder_eval, filterType_eval, bandwidth_eval, ARmethod, t_start, t_stop, targetPhase, maxRate, minNumPulses)

%%%% Phase-Locked Stimulation (Simulated) %%%%
%
% INPUTS
%   data - input data
%   f_s - sampling frequency
%   segmentLength - length of a signal segment, in sec
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   fstepsize - frequency step size for passband optimization
%   ARorder - order of autoregressive model
%   lambda - threshold fraction of power in optimized frequency band
%   filterOrder - order of bandpass filter
%   filterType - type of bandpass filter
%   ARorder_eval - AR order used to evaluate phase of simulated stimulation
%   filterOrder_eval - filter order used to evaluate phase of simulated 
%       stimulation
%   filterType_eval - filter type used to evaluate phase of simulated 
%       stimulation
%   bandwidth_eval - frequency bandwidth about central frequency of filter 
%       used to evaluate phase of simulated stimulation
%   ARmethod - method used to build AR model
%   t_start - start time of portion of signal used for AR model, in sec
%   t_stop - stop time of portion of signal used for AR model, in sec
%   targetPhase - (0 corresponds to peak, 180 corresponds to trough)
%   maxRate - rate of sliding window shifting, in Hz
%   minNumPulses - min number of output stimulation pulses required to
%       calculate phase-locking metrics
%
% OUTPUTS
%   indexStimTimes - index along 'data' of all simulation stimulation
%   stimFreqs - instataneous frequency at times of simulated stimulation
%   stimPhases - instataneous phase at times of simulated stimulation
%   phaseLockingMetrics - phase-locking metrics for all stimPhases
%       format is [mean phase, 95% CI lower of mean phase, 95% CI upper of
%       mean phase, circular variance, Rayleigh's test p-value, Rayleigh's
%       test z-score]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Example Parameter Initialization %%%%
% f_s = 500;
% segmentLength = 1; % 1 sec
% f_L = 4;
% f_H = 9;
% fstepsize = 0.05;
% ARorder_eval = 50;
% filterOrder_eval = 2;
% filterType_eval = 'butterworth';
% bandwidth_eval = 2;
% ARmethod = 'yw';
% targetPhase = 0;
% maxRate = 10;
% minNumPulses = 10;
% filters = {'butterworth', 'chebyshev1', 'chebyshev2', 'elliptic', 'bessel'};
% ARorder = optParams(i,1);
% lambda = optParams(i,2);
% filterOrder = optParams(i,3);
% filterType = filters{optParams(i,4)};
% t_stop = optParams(i,5);
% t_start = segmentLength - t_stop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

segmentSamples = f_s * segmentLength;
numSegments = floor(maxRate * (length(data)-2*segmentSamples) / f_s);

frequencies = f_L:fstepsize:f_H;

indexStimTimes = [];
stimFreqs = [];
stimPhases = [];
stimDelivered = 0;

for i=1:numSegments

    indexStart = ((i-1)*round(segmentSamples/maxRate)+1);
    indexStop = ((i-1)*round(segmentSamples/maxRate))+segmentSamples;
    signal_segment = data(indexStart:1:indexStop);
    Pxx = pyulear(signal_segment, ARorder, frequencies, f_s);
    
    if stimDelivered==1
        signal_eval = data((indexStimTimes(end)-f_s/2+1):1:(indexStimTimes(end)+f_s/2));        
        Pxx_eval = pyulear(signal_segment, ARorder_eval, frequencies, f_s);
        [~,f_cent_index_eval] = max(Pxx_eval);
        f_cent_eval = frequencies(f_cent_index_eval);
        if f_cent_eval-bandwidth_eval/2<f_L
            f_L_eval = f_L;
        else
            f_L_eval = f_cent_eval - bandwidth_eval/2;
        end
        if f_cent_eval+bandwidth_eval/2>f_H
            f_H_eval = f_H;
        else
            f_H_eval = f_cent_eval + bandwidth_eval/2;
        end
        [signal_eval_filtered] = zeroPhaseFilter(signal_eval, filterOrder_eval, filterType_eval, f_L_eval, f_H_eval, f_s);
        [phi_inst_eval, f_inst_eval] = instPhaseFreq(signal_eval_filtered, f_s/2, f_s);
        stimFreqs = [stimFreqs f_inst_eval];
        stimPhases = [stimPhases phi_inst_eval];
    end
    
    if ~isempty(findpeaks(Pxx)) && i<numSegments
        [fhat_L, fhat_H] = freqBandOpt(signal_segment, f_L, f_H, f_s, fstepsize, ARorder, lambda);
        [signal_segment_filtered] = zeroPhaseFilter(signal_segment, filterOrder, filterType, fhat_L, fhat_H, f_s);
        [signal_predicted, index0] = linearPredAR(signal_segment_filtered, ARorder, ARmethod, f_s, t_start, t_stop);
        [phi_inst, f_inst] = instPhaseFreq(signal_predicted, index0, f_s);
        t_delay = (1/f_inst) * mod(targetPhase-phi_inst+360,360) / 360;
        
        if t_delay>0 && t_delay<1/f_L        
            indexStimTimes = [indexStimTimes indexStop+round(t_delay*f_s)];
            stimDelivered = 1;
        else
            stimDelivered = 0;
        end
    else
        stimDelivered = 0;
    end

end

if length(indexStimTimes)>minNumPulses
    [mu ul ll] = circ_mean((stimPhases.*pi./180)');
    mu=mu.*180./pi;
    ul=ul.*180./pi;
    ll=ll.*180./pi;
    [var] = circ_var((stimPhases.*pi./180)');
    [pval, z] = circ_rtest((stimPhases.*pi./180)');
    phaseLockingMetrics = [mu,ll,ul,var,pval,z];
else
    phaseLockingMetrics = [nan,nan,nan,nan,nan,nan];
end

end


