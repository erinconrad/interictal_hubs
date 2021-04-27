function remove_empty_pts

%{
This is a possibly one time use function to remove extra entries from the
patient structure
%}

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
to_remove = zeros(length(pt),1);

% Loop over patients
for p =1:length(pt)
    
    name = pt(p).name;
    if isempty(name)
        to_remove(p) = 1;
    end
    
    
end

pt(logical(to_remove)) = [];
save([data_folder,'pt.mat'],'pt');

end