function plot_spikes_by_ch(whichPt,whichChs,whichFiles,userBlocks,spikes)

%% General parameters
do_car = 1; % common average ref? If 0, will do bipolar
surround = 5;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

% scripts and data
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];



%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
p = whichPt;

pt_name = pt(p).name;


%% Load spike file
if isempty(spikes)
    spike_folder = [results_folder,'new_spikes/'];
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
end


%% concatenate all spikes into one long thing
% Include an extra column for the file index
all_spikes = [];
all_elecs = {};

if isempty(whichFiles)
    whichFiles = 1:length(spikes.file);
end
for f = whichFiles
    if isempty(userBlocks)
        whichBlocks = 1:length(spikes.file(f).block);
    else
        whichBlocks = userBlocks;
    end
    for h = whichBlocks
        if isempty(spikes.file(f).block(h).gdf)
            continue;
        end
        all_spikes = [all_spikes;spikes.file(f).block(h).gdf,...
            repmat(f,size(spikes.file(f).block(h).gdf,1),1),...
            repmat(h,size(spikes.file(f).block(h).gdf,1),1)];
        all_elecs = [all_elecs;spikes.file(f).block(h).chLabels(spikes.file(f).block(h).gdf(:,1))];
    end
end

%% clean both sets of chlabels
whichChs = clean_labels_2(whichChs);
all_elecs = clean_labels_2(all_elecs);

%% Restrict spikes to those matching the electrodes of interest
matching_idx = ismember(all_elecs,whichChs);
all_spikes(~matching_idx,:) = [];
all_elecs(~matching_idx) = [];

%% Prep plot
which_plot = 0;
fig2 = figure;
set(gcf,'position',[0 0 1400 250])


while 1
    

    %% Randomly pick spike
    if isempty(all_spikes)
        fprintf('\nWarning, no spikes seen for this block and channel combo, choose another\n');
        close(fig2)
        break
    end
    sp = randi(size(all_spikes,1));

    %% Get info about the spike
    f = all_spikes(sp,3);
    h = all_spikes(sp,4);
    sp_time = all_spikes(sp,2);
    sp_ch = all_spikes(sp,1);
    tmul = spikes.file(f).block(1).params.tmul;
    absthresh = spikes.file(f).block(1).params.absthresh;
    fs = spikes.file(f).block(1).fs;
    fname = pt(p).ieeg.file(f).name;
    chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
    which_chs = spikes.file(f).block(h).run_chs;
    matching_chs = find(ismember(chLabels,whichChs));


    %% Get the EEG data
    %session = IEEGSession(fname, login_name, pwfile);
    run_times = [sp_time - surround,sp_time+surround];
    run_idx = run_times(1)*fs:run_times(2)*fs;
    %values = session.data.getvalues(run_idx,':');
    %session.delete;
    values = pull_ieeg_data(fname, login_name, pwfile, run_idx);
    sp_index = surround*fs;

    clean_labs = clean_labels_2(chLabels);
    if do_car
        values = new_pre_process(values,which_chs);
    else
        [values,bipolar_labels,chs_in_bipolar] = pre_process(values,clean_labs);
    end
    
        

    %% Plot data
    offset = 0;
    ch_offsets = zeros(length(matching_chs),1);
    ch_bl = zeros(length(matching_chs),1);
    for i = 1:length(matching_chs)
        ich = matching_chs(i);
        plot(linspace(0,surround*2,size(values,1)),values(:,ich)-offset,'linewidth',2);
        hold on
        ch_offsets(i) = offset;
        ch_bl(i) = -offset + nanmedian(values(:,ich));
        if ich == sp_ch
            sp_offset = offset;
        end
        if do_car
            text(surround*2+0.05,ch_bl(i),sprintf('%s',clean_labs{ich}),'fontsize',20)
        else
            text(surround*2+0.05,ch_bl(i),sprintf('%s',bipolar_labels{ich}),'fontsize',20)
        end
        if i<length(matching_chs)
            if ~isnan(min(values(:,ich)) - max(values(:,matching_chs(i+1))))
                offset = offset - (min(values(:,ich)) - max(values(:,matching_chs(i+1))));
            end
        end
    end
    
    plot(surround,values(round(sp_index),sp_ch)-sp_offset,'o','markersize',10)
    if do_car
        title(sprintf('Spike %d %1.1f s %s file %d block %d, tmul %d absthresh %d',...
            sp,sp_time,clean_labs{sp_ch},f,h,tmul,absthresh),'fontsize',10)
    else
        title(sprintf('Spike %d %1.1f s %s file %d block %d, tmul %d absthresh %d',...
            sp,sp_time,bipolar_labels{sp_ch},f,h,tmul,absthresh),'fontsize',10)
    end
    xlabel('Time (seconds)')
    set(gca,'fontsize',20)
    
    
    hold off

    str = input('\nPlot another? (y/n)\n','s');
    if contains(str,'y') || contains(str,'Y')
        % go to next
    else
        close(fig2)
        break % end function
    end

end
    


end