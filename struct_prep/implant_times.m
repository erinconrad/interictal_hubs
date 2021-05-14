%% Parameters
block = 60*30; % check every 30 minutes
mini_block = 60*30; % 30 minute block every thirty minutes (full block);

%% Be very careful about changing this
overwrite = 1; % overwrite if already exists? Would really screw up spike detections I already did

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'revision_pt.mat']);
pt = pt.pt;


whichPts = 1:length(pt);

for i = 1:length(whichPts)
    p = whichPts(i);
    name = pt(p).name;
    
    if isempty(pt(p).ieeg) || (isfield(pt(p).ieeg.file(1),'block') && overwrite == 0)
        fprintf('\nSkipping %s\n',name);
        continue;
    else
        fprintf('\nDoing %s\n',name);
    end
    
    
    
    if isempty(pt(p).ieeg), continue; end
    
    nfiles = length(pt(p).ieeg.file);
    for f = 1:nfiles
        
        % get duration
        dur = pt(p).ieeg.file(f).duration;

        % Split duration into chunks
        nblocks = ceil(dur/1e6/block); % 1e6 because duration is in microseconds
        bs = (0:block:floor(dur/1e6))';
        be = bs + block;
        be(end) = dur/1e6;
        pt(p).ieeg.file(f).block = [];
        pt(p).ieeg.file(f).block_times = [bs,be];
        for b = 1:nblocks
            pt(p).ieeg.file(f).block(b).start = bs(b);
            pt(p).ieeg.file(f).block(b).end = be(b);
            pt(p).ieeg.file(f).block(b).spikes = [];   
            pt(p).ieeg.file(f).block(b).run = [bs(b) be(b)];
        end


    end

end
    
save([data_folder,'revision_pt.mat'],'pt')
    
