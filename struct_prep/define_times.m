%% Parameters
block = 3600; % check every hour
mini_block = 60;

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


whichPts = 1:length(pt);

for i = 1:length(whichPts)
    p = whichPts(i);
    name = pt(p).name;
    
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

        pt(p).ieeg.file(f).block_times = [bs,be];
        for b = 1:nblocks
            pt(p).ieeg.file(f).block(b).start = bs(b);
            pt(p).ieeg.file(f).block(b).end = be(b);
            pt(p).ieeg.file(f).block(b).spikes = [];
            
            % pick a random minute in the hour long block
            s = randi([0 block-mini_block]);
            start_time = min(bs(b) + s,be(b) - mini_block);
            end_time = start_time + mini_block;
            
            run = [start_time,end_time];
            if run(1) < bs(b) || run(2) > be(b)
                if b == nblocks % may not have enough time to do the last bit
                    run = [];
                else
                    error('what');
                end
            end
            pt(p).ieeg.file(f).block(b).run = run;
        end


    end

end
    
save([data_folder,'pt.mat'],'pt')
    
