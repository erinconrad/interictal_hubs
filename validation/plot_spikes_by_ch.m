function plot_spikes_by_ch(whichPt,whichChs,do_car)

%% General parameters
surround = 5;

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
spike_folder = [results_folder,'new_spikes/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
p = whichPt;


pt_name = pt(p).name;


%% Load spike file
spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
spikes = spikes.spikes;

%% concatenate all spikes into one long thing
% Include an extra column for the file index
all_spikes = [];
all_elecs = {};
for f = 1:length(spikes.file)
    for h = 1:length(spikes.file(f).block)
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
figure
set(gcf,'position',[0 0 1400 250])


while 1

    %% Randomly pick spike
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
    chLabels = spikes.file(f).block(1).chLabels;
    which_chs = spikes.file(f).block(h).run_chs;



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
    plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
    hold on
    plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
    if do_car
        title(sprintf('Spike %d %1.1f s %s file %d, tmul %d absthresh %d',...
            sp,sp_time,clean_labs{sp_ch},f,tmul,absthresh),'fontsize',10)
    else
        title(sprintf('Spike %d %1.1f s %s file %d, tmul %d absthresh %d',...
            sp,sp_time,bipolar_labels{sp_ch},f,tmul,absthresh),'fontsize',10)
    end
    xlabel('Time (seconds)')
    set(gca,'fontsize',20)
    
    
    pause
    hold off


end
    


end