%{
To run any of the code, you will need to put a file in your path called
"interictal_hub_locations.m" that will point to various other paths. Here
is an example of what it should look like:
%}

%% Example locations file
%{
function locations = interictal_hub_locations

locations.script_folder = [locations.main_folder,'scripts/']; (path to this
code base)
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




%}