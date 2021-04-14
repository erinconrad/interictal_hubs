function plot_example_detections(whichPts)

%% General parameters
n_sp = 50;
n_per_fig = 10;
surround = 5;

%% Locations
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
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
    out_folder = [results_folder,'validation/',pt_name,'/'];
    if exist(out_folder,'dir') == 0, mkdir(out_folder); end

    %% Load spike file
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    
    %% concatenate all spikes into one long thing
    % Include an extra column for the file index
    all_spikes = [];
    for f = 1:length(spikes.file)
        for h = 1:length(spikes.file(f).hour)
            all_spikes = [all_spikes;spikes.file(f).hour(h).gdf,...
                repmat(f,size(spikes.file(f).hour(h).gdf,1),1)];
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
        sp_label = spikes.file(f).hour(1).chLabels{sp_ch};
        tmul = spikes.file(f).hour(1).params.tmul;
        absthresh = spikes.file(f).hour(1).params.absthresh;
        fs = spikes.file(f).hour(1).fs;
        fname = pt(p).ieeg.file(f).name;
        
        %% Get the EEG data
        session = IEEGSession(fname, login_name, pwfile);
        run_times = [sp_time - surround,sp_time+surround];
        run_idx = run_times(1)*fs:run_times(2)*fs;
        values = session.data.getvalues(run_idx,':');
        sp_index = surround*fs;
        
        %% Plot data
        axes(ha(b))
        plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
        hold on
        plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
        title(sprintf('Spike %d %1.1f s %s file %d, tmul %d absthresh %d',...
            sp,sp_time,sp_label,f,tmul,absthresh),'fontsize',10)
        if b ~= n_per_fig
            xticklabels([])
        end
        yticklabels([])
        set(gca,'fontsize',10)

        if b == n_per_fig
            xlabel('Time (s)')
            which_plot = which_plot + 1;
            print([out_folder,sprintf('spikes_%d',which_plot)],'-dpng');
            close(gcf)
        end
        
    end
    
end

end