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