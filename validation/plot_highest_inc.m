function plot_highest_inc(out)

% note there is currently a line of code that will break this if I try to
% do lowest channels instead

surround = 24;
nhighest = 10;
nperch = 5;
sp_sur = 1.5;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'rate_validation/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

for p = 1:length(out)
    
    %% Prep figure
    figure
    set(gcf,'position',[100 100 1200 1000])
    tiledlayout(nhighest,nperch,'padding','tight','tilespacing','tight');
    
    
    %% Get rate info
    name = out(p).name;
    
    pt_p = nan;
    for j = 1:length(pt)
        if strcmp(pt(j).name,name)
            pt_p = j;
            break
        end
    end


    % load spikes
    spikes = load([spike_folder,sprintf('%s_spikes.mat',name)]);
    spikes = spikes.spikes;
    
    
    rate = out(p).rate;
    chLabels = out(p).unchanged_labels;
    cblock = out(p).change_block;

    ekg = identify_ekg_scalp(chLabels);

    rate = rate(~ekg,:);
    chLabels = chLabels(~ekg);
    findices = out(p).findices;
    bindices = out(p).bindices;

    [pre,post] = get_surround_times(rate,cblock,surround);
    
    
    % Get bindices and findices corresponding to post
    post_b = bindices(post);
    post_f = findices(post);
    
    if length(unique(post_f)) ~= 1
        error('why');
    end
    post_f = post_f(1);
    
    % Get relative rate change
    pre_rate = nanmean(rate(:,pre),2);
    post_rate = nanmean(rate(:,post),2);
    rel_rate_change = (post_rate-pre_rate)./abs(pre_rate);
    
    rel_rate_change(isnan(rel_rate_change)) = -inf;
    [~,I] = sort(rel_rate_change,'descend');
    
    % Take n highest
    highest_idx = I(1:nhighest);
    highest_labels = chLabels(highest_idx);
    
    %% Get the appropriate spikes
    % concatenate all spikes into one long thing
    % Include an extra column for the file index
    all_spikes = [];
    all_elecs = {};

    for h = post_b
        if isempty(spikes.file(post_f).block(h).gdf)
            continue;
        end
        
        all_spikes = [all_spikes;spikes.file(post_f).block(h).gdf,...
            repmat(post_f,size(spikes.file(post_f).block(h).gdf,1),1),...
            repmat(h,size(spikes.file(post_f).block(h).gdf,1),1)];
        all_elecs = [all_elecs;spikes.file(post_f).block(h).chLabels(spikes.file(post_f).block(h).gdf(:,1))];
    end
    
    % clean both sets of chlabels
    whichChs = clean_labels_2(highest_labels);
    all_elecs = clean_labels_2(all_elecs);

    for ich = 1:nhighest
        curr_label = whichChs{ich};
        matching_idx = ismember(all_elecs,curr_label);
        
        curr_spikes = all_spikes(matching_idx,:);
        curr_elecs = all_elecs(matching_idx);
        nspikes = size(curr_spikes,1);
        
        num_to_plot = min([nspikes,nperch]);
        
        % Choose random spikes
        rand_spikes = randsample(nspikes,num_to_plot);
        
        %% Get spike
        for is = 1:num_to_plot
            
            % get the spike index
            sp = rand_spikes(is);
            
            % Info about the spike
            f = curr_spikes(sp,3);
            h = curr_spikes(sp,4);
            sp_time = curr_spikes(sp,2);
            fs = spikes.file(f).block(1).fs;
            fname = pt(pt_p).ieeg.file(f).name;
            which_chs = spikes.file(f).block(h).run_chs;
            chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
            sp_index = sp_sur*fs;
            
            fprintf('\nShowing spike %d on %s at %1.1fs, file %d, block %d\n',...
                sp,curr_label,sp_time,f,h);
            
            % Re-derive appropriate sp_ch
            sp_ch = find(ismember(chLabels,curr_label));
            
            % Get the EEG data
            run_times = [sp_time - sp_sur,sp_time+sp_sur];
            run_idx = run_times(1)*fs:run_times(2)*fs;
            values = pull_ieeg_data(fname, login_name, pwfile, run_idx);
            clean_labs = clean_labels_2(chLabels);
            values = new_pre_process(values,which_chs);
            
            % Plot data
            nexttile
            plot(linspace(0,sp_sur*2,size(values,1)),values(:,sp_ch),'linewidth',2);
            hold on
            plot(sp_sur,values(round(sp_index),sp_ch),'o','markersize',10)
            %{
            title(sprintf('Spike %d %1.1f s %s file %d block %d',...
            sp,sp_time,clean_labs{sp_ch},f,h),'fontsize',10)
            %}
            %xlabel('Time (seconds)')
            %set(gca,'fontsize',20)
            xticklabels([])
            yticklabels([])
            if is == 1
                ylabel(clean_labs{sp_ch});
            end
            
        end
        
        % Pad few spike cases
        for k = 1:nperch-num_to_plot
            nexttile
        end
        
    end
    
    % Save file
    print(gcf,[out_folder,name],'-dpng')
    close(gcf)
    
end

end