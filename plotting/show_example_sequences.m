function show_example_sequences(whichPts)

%% General parameters
do_simple = 1;
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


if isempty(whichPts)
    listing = dir([spike_folder,'*.mat']);
    for i = 1:length(listing)
        C = listing(i).name;
        temp_name = strsplit(C,'_');
        temp_name = temp_name{1};
        for j = 1:length(pt)
            pt_name = pt(j).name;
            if strcmp(temp_name,pt_name)
                whichPts = [whichPts,j];
                break
            end
        end
    end
end

for p = whichPts
    
    pt_name = pt(p).name;
    
     %% Load spike file
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    
    nfiles = length(spikes.file);
    for f = 1:nfiles
        
        fname = pt(p).ieeg.file(f).name;
        chLabels = spikes.file(f).block(1).chLabels;
        nchs = length(spikes.file(f).block(1).chLabels);
        old_gdf = [];
        %% Get old gdf
        for b = 1:length(spikes.file(f).block)
            old_gdf = [old_gdf;spikes.file(f).block(b).gdf];
        end
        [~,~,old_seq] = get_sequences(old_gdf,nchs);
        
        %% Get spike details as a table
        T = convert_details_to_table(spikes,f,[]);
        
        %% Get sequences
        times = T.peak_idx./T.fs + T.run_start;
        [times,I] = sort(times);
        T = T(I,:);
        chs = T.ch;
        gdf = [chs,times];
        seq = new_get_sequences(gdf);
        
        nseq = length(seq);
        
        
        while 1
            %% Choose a random sequence
            s = randi(nseq);
            curr_seq = seq{s};
            which_spikes = curr_seq(:,3);
            curr_indices = T.peak_idx(which_spikes);
            rise = T.mid_rise_idx(which_spikes);
            fall = T.mid_fall_idx(which_spikes);
            curr_chs = T.ch(which_spikes);
            fs = T.fs(which_spikes(1));
            run_start = T.run_start(which_spikes(1));

            %% Get run times
            start_time = run_start + curr_indices(1)/fs - surround;
            end_time = run_start + curr_indices(1)/fs + surround;
            run_times = [start_time end_time];
            dur = surround*2;

            %% Get ieeg
            run_idx = round(run_times(1)*fs):round(run_times(2)*fs);
            values = pull_ieeg_data(fname, login_name, pwfile, run_idx);


            %% Adjust indices to this new start time
            indices = curr_indices - curr_indices(1) + surround*fs;
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

            for i = 1:length(indices)
                %index = spikes(s,1);

                % convert index to time
                time = indices(i)/fs;
                rtime = rise(i)/fs;
                ftime = fall(i)/fs;

                ch = curr_chs(i);
                offset_sp = ch_offsets(i);

                plot(time,values(indices(i),ch) - offset_sp,'bo')
                plot(rtime,values(rise(i),ch) - offset_sp,'go')
                plot(ftime,values(fall(i),ch) - offset_sp,'ro')

            end
            
            title(sprintf('%s file %d time %1.1f',pt_name,f,run_times(1)))
            pause
            close(gcf)
        end

    
    
end


end