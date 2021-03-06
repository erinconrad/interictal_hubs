%% Implant Effect ReadMe

%{
This folder contains the scripts to run the analysis for the implant effect
paper. Note that some scripts that are shared with other project exist
outside of this folder in the main interictal hubs codebase, and so the
entire interictal hubs code base should be downloaded in order to run this
analysis. 


To replicate the Implant Effect analyses, do the following:

1) Download the interictal_hubs codebase, located at:
https://github.com/erinconrad/interictal_hubs

2) Add a file to your Matlab path titled "interictal_hub_locations.m". This
points to various other paths.

Here is what this file should look like:

"function locations = interictal_hub_locations

locations.main_folder = [path to where you will output results from the
analysis]
locations.script_folder = [path to this code base]
locations.ieeg_folder = [path to the ieeg.org codebase]
locations.ieeg_pw_file = [path to you ieeg.org password file]
locations.ieeg_login = [your ieeg.org login name]

end"

You can leave the ieeg information empty if you will just be running the
analyses (whereas you will need the ieeg code base if you want to
re-calculate spikes, Pearson connectivity networks, or alpha-delta ratios,
as described at the bottom of this file).

3) In the scripts/implant_analyses/ folder (the folder containing this
ReadMe), run run_all_analysis.m:

>> run_all_analysis

This will call functions in the scripts/implant_analyses/main_analyses/
folder to run the different analyses from the paper.

4) If you wish to re-calculate intermediate data (spikes, Pearson networks,
or alpha-delta ratios) you will need access to the IEEG codebase, because
this involves downloading raw EEG data from ieeg.org. You will first need
to download the IEEG toolbox at https://www.ieeg.org/ and create a free
account.

The path to the code to calculate spikes is get_hubs/new_spikes.m
The path to the code to calculate PC networks is get_hubs/pc_networks.m
The path to the code to calculate alpha/delta ratios is get_hubs/alpha_delta_ratio.m


Erin Conrad
University of Pennsylvania
July 2021


%}