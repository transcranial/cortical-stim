function [signal_predicted, index0] = linearPredAR(signal, ARorder, method, f_s, t_start, t_stop)

%%%% Autoregressive Time-Series Forward Prediction %%%%
%
% Builds AR model of ARORDER by METHOD and forecasts (double) time series
% SIGNAL N steps into the future. Returns forecasted (double) time series YF of
% length K.
%
% INPUTS
%   signal - input signal used as basis for building AR model
%   ARorder - order of AR model
%   method - method for building AR model
%   f_s - sampling frequency
%   t_start - start time of portion of signal used for AR model, in sec
%   t_stop - stop time of portion of signal used for AR model, in sec
%
% OUTPUTS
%   signal_predicted - AR predicted signal that is twice the length between
%       the point corresponding to t_stop and t=0 (t=0 is at the end of the
%       the original signal)
%   index0 - index of signal_predicted corresponding to t=0
%
% METHOD can be:
% 1) 'burg' (Burg's lattice-based method solving the lattice filter
% equations using the harmonic mean of forward and backward squared
% prediction errors),
% 2) 'fb' (Forward-backward approach minimizing the sum of a least-squares
% criterion for a forward model, and the analogous criterion for a time-
% reversed model),
% 3) 'gl' (Geometric lattice approach similar to Burg's method but uses the
% geometric mean instead of the harmonic mean during minimization),
% 4) 'ls' (Least-squares approach minimizing the standard sum of squared
% forward-prediction errors),
% 5) 'yw' (Yule-Walker approach solving the Yule-Walker equations, formed
% from sample covariances)
%
% For t_start and t_stop, as an example, if signal is 1 sec in length, then
% t_start=0.25, t_stop=0.75 uses the middle 0.5 sec of signal to build AR
% model. The predicted signal will then be 0.5 second in length,
% corresponding to the time from t_stop (t=-0.25s) to t=0.25s.
%
%%%%

k = 2 * (length(signal) - round(f_s*t_stop));
index0 = k/2;

% Process input data
y = iddata(signal(round(f_s*t_start):round(f_s*t_stop))');
[N, ny, nu] = size(y);
zer = iddata(zeros(k,ny),zeros(k,nu),y.Ts);
yz = [y; zer];

% Build model
mb = ar(y,ARorder,method);

% Perform prediction
yft = predict(mb,yz,k,'e');

% Keep only predicted time series
yff = yft(N+1:N+k);
signal_predicted = yff.OutputData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end