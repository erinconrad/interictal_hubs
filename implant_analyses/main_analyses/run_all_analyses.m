%{
This script calls the functions to run all major analyses for the implant
paper. It assumes that you have an file called "out.mat" containing spike
times, Pearson correlation networks, and alpha delta ratios for each
patient.
%}


%% Set random generator seed
% This is so that I get the same result (assuming I don't change the code)
% each time I run this. Because I am doing a Monte Carlo test there is some
% randomness.
rng(100)

%% Run spike and node strength stability analysis
new_ros_fig([],1);

%% Run distance-rate change correlation analysis
all_corrs([],1);

%% Run distance-FC-cosi correlation analysis
dist_cosi_pc([],1);

%% Run rate analysis
rate_all_pts([],1);

%% Run alpha-delta ratio analysis
ad_analyses([],1);

%% Run anatomy analysis
anatomy_analyses([],1);