%% Interictal hubs ReadMe
%{
To run analyses for the Implant Effect paper, navigate to the
implant_analysis/ folder and read the ReadMe in that folder.


%}


%{
Notes:
- results/backups/spikes_revision_August17/ contains the _spikes.mat files
used for the implant revision project as of August 17 (pre implementing new
detector)

I made all sorts of crazy stuff to fix the jackbox switches for Hup132. To
summarize, I made code at the get_hubs stage to redefine channel labels. I
then made code to break the output into new files so that ch labels are
consistent within a file. I then edited the get_gdf_details code to take
locs and anatomy from the spikes file rather than the pt file (because the
pt file still retains the original file definitions and so is wrong!). 

I still need to fix this for the pc and ad stuff and I need to check that
the analysis files don't call the pt structure (which would be wrong for
this patient).

%}

%{
To run any of the code, you will need to put a file in your path called
"interictal_hub_locations.m" that will point to various other paths. Here
is an example of what it should look like:
%}

%% Example locations file
%{
function locations = interictal_hub_locations

locations.main_folder = [path to where you will output results from the
analysis]
locations.script_folder = [path to this code base]
locations.ieeg_folder = [path to the ieeg.org codebase]
locations.ieeg_pw_file = [path to you ieeg.org password file]
locations.ieeg_login = [your ieeg.org login name]

end
%}

%{ 
You will also need the ieeg.org code base and you need to point to it in
the file described above
%}

%% Figure and analyses scripts for implant effect
%{
get_gdf_details: this function takes the Name_spikes.mat files and gets
details like spike rate over time and co-spikiness, along with information
like distance from newest added electrodes, locations, sampling rate, index
of the block with the revision, etc. It saves an out.mat structure with
this information

rate_all_pts: this generates Figure 1 for the implant effect paper, showing
the spike rate over time for a single patient and showing the change in
spike rate pre- to post-revision for all patients.

ros_all_pts: this generates Figure 2 for the implant effect paper, showing
a raster plot of spike rate and ns for a single patient and showing the
rate order stability and NS order stability from pre- to post-revision for
all patients.

corr_each_pt: this generates Figure 3 for the implant paper, showing a
scatter plot of rate change as a function of something for each patient.


%}