function [optParams, indexStimTimes, stimFreqs, stimPhases, phaseLockingMetrics] = findOptParams(data, f_s)
warning off;

%%%% Offline Genetic Algorithm for Optimizing Algorithm Parameters %%%%
%
% INPUTS
%   data - input data
%   f_s - sampling frequency
%
% OUTPUTS
%   optParams - optimized parametres in format of: [ARorder, lambda,
%       filterOrder, filterType, t_stop]. t_start is segmentLength - t_stop.
%   indexStimTimes - index along 'data' of all simulation stimulation
%   stimFreqs - instataneous frequency at times of simulated stimulation
%   stimPhases - instataneous phase at times of simulated stimulation
%   phaseLockingMetrics - phase-locking metrics for all stimPhases
%       format is [mean phase, 95% CI lower of mean phase, 95% CI upper of
%       mean phase, circular variance, Rayleigh's test p-value, Rayleigh's
%       test z-score]
%%%%

% Preset parameters:
%   segmentLength - length of a signal segment, in sec
%   f_L - passband lower frequency limit
%   f_H - passband higher frequency limit
%   fstepsize - frequency step size for passband optimization
%   ARorder_eval - AR order used to evaluate phase of simulated stimulation
%   filterOrder_eval - filter order used to evaluate phase of simulated 
%       stimulation
%   filterType_eval - filter type used to evaluate phase of simulated 
%       stimulation
%   bandwidth_eval - frequency bandwidth about central frequency of filter 
%       used to evaluate phase of simulated stimulation
%   ARmethod - method used to build AR model
%   targetPhase - (0 corresponds to peak, 180 corresponds to trough)
%   maxRate - rate of sliding window shifting, in Hz
%   minNumPulses - min number of output stimulation pulses required to
%       calculate phase-locking metrics
segmentLength = 1; % 1 sec
f_L = 4;
f_H = 9;
fstepsize = 0.05;
ARorder_eval = 50;
filterOrder_eval = 2;
filterType_eval = 'butterworth';
bandwidth_eval = 2;
ARmethod = 'yw';
targetPhase = 0;
maxRate = 10;
minNumPulses = 10;

% Undetermined parameters: 
% 	ARorder
%   lambda
%   filterOrder
%   filterType
%   t_start, t_stop

filters = {'butterworth', 'chebyshev1', 'chebyshev2', 'elliptic', 'bessel'};
lb = [2;0.5;0.01;0.01;0.55];
ub = [100;1;5;5;0.95];
options = gaoptimset('Display','iter','CreationFcn',@gacreationlinearfeasible,...
    'PopInitRange',[lb,ub]','Generations',50,'PopulationSize',200,'EliteCount',20,...
    'CrossoverFraction',0.7,...
    'FitnessScalingFcn',@fitscalingrank,'MutationFcn',@mutationadaptfeasible,...
    'SelectionFcn',@selectionroulette,'CrossoverFcn',@crossoverheuristic,...
    'StallGenLimit',5, 'TolFun',1e-4);
[x,~,~,~] = ga(@fitnessfcn,5,[],[],[],[],lb,ub,[],options);

optParams = [round(x(1)), round(x(2)*100)/100, ceil(x(3)), ceil(x(4)), round(x(5)*100)/100];
ARorder = round(x(1));
lambda = round(x(2)*100)/100;
filterOrder = ceil(x(3));
filterType = filters{ceil(x(4))};
t_stop = round(x(5)*100)/100;
t_start = segmentLength - t_stop;
[indexStimTimes, stimFreqs, stimPhases, phaseLockingMetrics] = stimSim(data, f_s, ...
    segmentLength, f_L, f_H, fstepsize, ARorder, lambda, filterOrder, filterType, ...
    ARorder_eval, filterOrder_eval, filterType_eval, bandwidth_eval, ARmethod, ...
    t_start, t_stop, targetPhase, maxRate, minNumPulses);

    function y = fitnessfcn(x)
        
        ARorder = round(x(1));
        lambda = round(x(2)*100)/100;
        filterOrder = ceil(x(3));
        filterType = filters{ceil(x(4))};
        t_stop = round(x(5)*100)/100;
        t_start = segmentLength - t_stop;
        
        [~,~,~,phaseLockingMetrics_temp] = stimSim(data, f_s, segmentLength, ...
            f_L, f_H, fstepsize, ARorder, lambda, filterOrder, filterType, ...
            ARorder_eval, filterOrder_eval, filterType_eval, bandwidth_eval, ...
            ARmethod, t_start, t_stop, targetPhase, maxRate, minNumPulses);
        y = abs(targetPhase - phaseLockingMetrics_temp(:,1))./180 + abs(phaseLockingMetrics_temp(:,3)-phaseLockingMetrics_temp(:,2))./180 + ARorder/100;
        if phaseLockingMetrics_temp(:,4)==0
            y=y+1;
        end
    end

end