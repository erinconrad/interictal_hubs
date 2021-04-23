function new_example_sequences(p)

%% General parameters
filt = 2;
which_files = [];
timing = 'peak';
surround = 2;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
addpath(genpath(locations.ieeg_folder));
spike_folder = [results_folder,'spikes/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


pt_name = pt(p).name;

 %% Load spike file
spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
spikes = spikes.spikes;

nfiles = length(spikes.file);
if isempty(which_files)
    which_files = 1:nfiles;
end

all_seq = {};

% Loop over files
for i = 1:length(which_files)
    f = which_files(i);

    fname = pt(p).ieeg.file(f).name;
    chLabels = spikes.file(f).block(1).chLabels;
    nchs = length(chLabels);
    nblocks = length(spikes.file(f).block);

    % Loop over blocks
    for h = 1:nblocks
        block = spikes.file(f).block(h);

        if block.run_skip == 1 || isempty(block.gdf)
            continue;
        end

        %% Get spike channels and times (whatever time you want)
        chs = block.details.filter(filt).gdf(:,1);
        times = block.details.filter(filt).(timing);
        peak = block.details.filter(filt).peak;
        rise = block.details.filter(filt).rise;
        fall = block.details.filter(filt).fall;
        fs = block.fs;

        %% Get sorted spike indices and chs
        [times,I] = sort(times);
        chs = chs(I);

        %% Construct gdf
        gdf = [chs,times,...
            repmat(f,length(chs),1),...
            repmat(h,length(chs),1),...
            peak,rise,fall];

        %% Get sequences
        seq = new_get_sequences(gdf,nchs,fs);

        %% put in all sequence cell
        all_seq = [all_seq;seq]; 

    end

end

nseq = length(all_seq);

while 1
    %% Choose a random sequence
    s = randi(nseq);
    curr_file = all_seq{s}(:,3);
    curr_block = all_seq{s}(:,4);
    curr_indices = all_seq{s}(:,2);
    curr_chs = all_seq{s}(:,1);
    peak = all_seq{s}(:,5);
    rise = all_seq{s}(:,6);
    fall = all_seq{s}(:,7);
    
    chLabels = spikes.file(curr_file(1)).block(1).chLabels;
    fs = spikes.file(curr_file(1)).block(1).fs;
    run_start = spikes.file(curr_file(1)).block(curr_block(1)).run_times(1);
    

    %% Get run times
    start_time = run_start + curr_indices(1)/fs - surround;
    end_time = run_start + curr_indices(1)/fs + surround;
    run_times = [start_time end_time];
    dur = surround*2;

    %% Get ieeg
    run_idx = round(run_times(1)*fs):round(run_times(2)*fs);
    values = pull_ieeg_data(fname, login_name, pwfile, run_idx);


    %% Adjust indices to this new start time
    peak = peak - curr_indices(1) + surround*fs;
    rise = rise - curr_indices(1) + surround*fs;
    fall = fall - curr_indices(1) + surround*fs;

    %% Plot
    figure
    offset = 0;
    ch_offsets = zeros(length(curr_chs),1);
    ch_bl = zeros(length(curr_chs),1);

    for i = 1:length(curr_chs)
        ich = curr_chs(i);
        plot(linspace(0,dur,size(values,1)),values(:,ich)-offset)
        ch_offsets(i) = offset;
        ch_bl(i) = -offset + nanmedian(values(:,ich));
        hold on
        text(dur+0.05,ch_bl(i),sprintf('%s',chLabels{ich}))

        if i<length(curr_chs)
            if ~isnan(min(values(:,ich)) - max(values(:,curr_chs(i+1))))
                offset = offset - (min(values(:,ich)) - max(values(:,curr_chs(i+1))));

            end
        end

    end

    for i = 1:length(peak)

        % convert index to time
        ptime = peak(i)/fs;
        rtime = rise(i)/fs;
        ftime = fall(i)/fs;

        ch = curr_chs(i);
        offset_sp = ch_offsets(i);

        plot(ptime,values(peak(i),ch) - offset_sp,'bo')
        plot(rtime,values(rise(i),ch) - offset_sp,'go')
        plot(ftime,values(fall(i),ch) - offset_sp,'ro')

    end

    title(sprintf('%s file %d time %1.1f',pt_name,f,run_times(1)))
    pause
    close(gcf)
end

    
    



end