# Method for Real-Time Cortical Oscillation Detection and Phase-Locked Stimulation

MATLAB source files for:

[Link](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5705563) Chen, L. Leon, et al. "Real-time brain oscillation detection and phase-locked stimulation using autoregressive spectral estimation and time-series forward prediction." Biomedical Engineering, IEEE Transactions on 60.3 (2013): 753-762.

[Link](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6090843) Chen, L. Leon, et al. "A method for real-time cortical oscillation detection and phase-locked stimulation." Engineering in Medicine and Biology Society, EMBC, 2011 Annual International Conference of the IEEE. IEEE, 2011.


### main functions

##### findOptParams 
offline genetic algorithm for optimizing algorithm parameters

##### findBestChannel 
finds best single channel based on: power in frequency band, temporal coherence of bandpass filtered signal, and by both metrics combined

##### stimSim 
simulated stimulation over an entire dataset (output is summary of stimulation times, frequency, phase, as well as summary of phase-locking metrics)

##### stimRT 
real-time stimulation (output is `t_delay`)

### auxillary support functions

##### freqBandOpt

##### instPhaseFreq

##### linearPredAR

##### zeroPhaseFilter

##### temporalCoherence