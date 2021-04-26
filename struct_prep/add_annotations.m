function add_annotations

%% Parameters
overwrite = 1; % overwrite if already exists?

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


whichPts = 1:length(pt);


for p = 1:length(pt)

for f = 1:length(pt(p).ieeg.file)
    
    clear event
    
     %% Download data
    loginname = 'erinconr';
    session = IEEGSession(pt(p).ieeg_names{f}, loginname, pwname);    
    
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
        pt(p).filename(f).ann(ai) = ann;
    end
    
    session.delete;
    
end


end

%% Save the file
save([data_folder,'pt.mat'],'pt');

end