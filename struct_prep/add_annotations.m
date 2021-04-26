function add_annotations

%% Parameters
overwrite = 0; % overwrite if already exists?

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(ieeg_folder))

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


whichPts = 1:length(pt);


for p = 1:length(pt)

for f = 1:length(pt(p).ieeg.file)
    
    if overwrite == 0
        if isfield(pt(p).ieeg.file(f),'ann') && ~isempty(pt(p).ieeg.file(f).ann)
            fprintf('\nSkipping %s file %d\n',pt(p).name,f);
            continue
        end
    end
            
    clear event
    
     %% Download data
    session = IEEGSession(pt(p).ieeg.file(f).name, login_name, pwfile);    
    
    %if f == 4, error('look\n'); end
    
    %% Get annotations
    n_layers = length(session.data.annLayer);
    for ai = 1:n_layers
        a=session.data.annLayer(ai).getEvents(0);
        n_ann = length(a);
        for i = 1:n_ann
            event(i).start = a(i).start/(1e6);
            event(i).stop = a(i).stop/(1e6); % convert from microseconds
            event(i).type = a(i).type;
            event(i).description = a(i).description;
        end
        ann.event = event;
        ann.name = session.data.annLayer(ai).name;
        pt(p).ieeg.file(f).ann(ai) = ann;
    end
    
    session.delete;
    
end


end

%% Save the file
save([data_folder,'pt.mat'],'pt');

end