% This is the main file,  which
% 1) allows the user to specify the main % parameters and the statistical model 
% 2) runs the simulations by calling the routine run_tests and save the results
% 3) creates Figure 1 and 2
% This is the only file which should be modifed by the user

close all
clear
clc
rng(1) % set the random number generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parameters that can be chosen by the user %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_name = 'scratch'; % the results will be saved in a file 'save_name.mat'

%%% Sample sizes to be tested (m denotes the number of signals in the training sample)
Nmm = 2;        % how many different values of m we try 
m_min = 3000;   % minimum sample size
m_max = 300000; % maximum sample size

%%% Experiment repetitions (to verify the generalization bound, which are in 
% expectation, we replicate them K times and take the average)
K = 30;

%%% Signal resolution: all the signals are supposed to be in L^2([0,1]) but
% are represented by dividing the interval in N pixels
NumbN = 2;  % how many different resolutions we test
Nmin = 64;  % minimum resolution (if using wavelets, must be a power of 2)
Nmax = 256; % maximum resolution (if using wavelets, must be a power of 2)
Nk = 1;     % resolution to be depicted in Figure 1 (a number between 1 and NumbN)

%%% Statistical model
% x is a Gaussian random variable with mean mu and
% covariance operator described by a convolution operator whose kernel k_Sx 
% is supported in (-c_k,c_k) being c_k < 1
mu = @(x) 1-abs(2*x-1);  % must be a handle-function
c_k = 0.2;               % kernel support size
k_Sx   = @(x) (1-exp(-abs(c_k./x).^4)); % must be a handle-function

% epsilon is either a Gaussian white noise or a uniform white noise, either
% wrt to the pixel basis or the wavelet components. Its mean is 0 and its
% covariance is sigma*I
sigma = 0.05;
error_model = 1; % 1: gaussian white noise
                 % 2: uniform white noise
                 % 3: uniform white noise wrt wavelet components
 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% No user's modifications from here on %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
run_tests; % script which contains the main routine

saved_plot_fig1(save_name,Nk) % creates Figure 1
saved_plot_fig2(save_name)    % creates Figure 2

                 
                 
               
