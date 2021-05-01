function plot_example_detections(whichPts,which_ver,overwrite)

%% General parameters
n_sp = 50;
n_per_fig = 10;
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
if which_ver == 1
    spike_folder = [results_folder,'spikes/'];
elseif which_ver == 2
    spike_folder = [results_folder,'new_spikes/'];
end

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
    out_folder = [results_folder,'validation/',pt_name,'/'];
    
    if exist(out_folder,'dir') ~= 0
        listing = dir([out_folder,'*.jpg']);
        if length(listing) == 5
            if overwrite == 0
                fprintf('\nSkipping %s\n',pt_name);
                continue
            end
        else
            fprintf('\nDoing %s\n',pt_name);
        end
    else
        fprintf('\nDoing %s\n',pt_name);
        mkdir(out_folder);
    end

    %% Load spike file
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    
    %% concatenate all spikes into one long thing
    % Include an extra column for the file index
    all_spikes = [];
    for f = 1:length(spikes.file)
        for h = 1:length(spikes.file(f).block)
            all_spikes = [all_spikes;spikes.file(f).block(h).gdf,...
                repmat(f,size(spikes.file(f).block(h).gdf,1),1)];
        end
    end
    
    which_plot = 0;
    for i = 1:n_sp

        b = mod(i,n_per_fig);
        if b == 1
            figure
            set(gcf,'position',[0 0 1400 800])
            [ha,~] = tight_subplot(n_per_fig,1,[0.02 0.02],[0.05 0.05],[0.02 0.02]);
        elseif b == 0
            b = n_per_fig; 
        end
        
        %% Randomly pick spike
        sp = randi(size(all_spikes,1));
        
        %% Get info about the spike
        f = all_spikes(sp,3);
        sp_time = all_spikes(sp,2);
        sp_ch = all_spikes(sp,1);
        tmul = spikes.file(f).block(1).params.tmul;
        absthresh = spikes.file(f).block(1).params.absthresh;
        fs = spikes.file(f).block(1).fs;
        fname = pt(p).ieeg.file(f).name;
        chLabels = spikes.file(f).block(1).chLabels;
        
        
        
        %% Get the EEG data
        %session = IEEGSession(fname, login_name, pwfile);
        run_times = [sp_time - surround,sp_time+surround];
        run_idx = run_times(1)*fs:run_times(2)*fs;
        %values = session.data.getvalues(run_idx,':');
        %session.delete;
        values = pull_ieeg_data(fname, login_name, pwfile, run_idx);
        sp_index = surround*fs;
        
        clean_labs = clean_labels_2(chLabels);
        [values,bipolar_labels] = pre_process(values,clean_labs);
        
        %% Plot data
        axes(ha(b))
        plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
        hold on
        plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
        title(sprintf('Spike %d %1.1f s %s file %d, tmul %d absthresh %d',...
            sp,sp_time,bipolar_labels{sp_ch},f,tmul,absthresh),'fontsize',10)
        if b ~= n_per_fig
            xticklabels([])
        end
        yticklabels([])
        set(gca,'fontsize',10)

        if b == n_per_fig
            xlabel('Time (s)')
            which_plot = which_plot + 1;
            print([out_folder,sprintf('spikes_%d',which_plot)],'-djpeg');
            close(gcf)
        end
        
    end
    
end

end