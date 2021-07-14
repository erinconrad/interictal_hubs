%{
This script calls the functions to run all major analyses for the implant
paper. It assumes that you have an file called "out.mat" containing spike
times, Pearson correlation networks, and alpha delta ratios for each
patient. This file should be in the implant_analyses/ folder (the same
folder as this script).
%}


%% Set random generator seed
% This is so that I get the same result (assuming I don't change the code)
% each time I run this. Because I am doing a Monte Carlo test there is some
% randomness.
rng(100)

%% Temporarily add this codebase to the Matlab path
locations = interictal_hub_locations; % This is a script you will need to add to your own path (see ReadMe for details)
addpath(genpath(locations.script_folder)); % path to the interictal_hubs codebase

%% Load the out file
% This file contains the timing of spike detections, the calculated Pearson
% correlation networks, the alpha delta ratios, and other assorted info
% (electrodes labels, electrode locations, anatomy). This should be in the
% same folder as this script (run_all_analyses)
out = load('out.mat');
out = out.out;

%% Run spike and node strength stability analysis
new_ros_fig([],1,out);

%% Run distance-rate change correlation analysis
all_corrs([],1,out);

%% Run distance-FC-cosi correlation analysis
dist_cosi_pc([],1,out);

%% Run rate analysis
rate_all_pts([],1,out);

%% Run alpha-delta ratio analysis
ad_analyses([],1,out);

%% Run anatomy analysis
anatomy_analyses([],1,out);