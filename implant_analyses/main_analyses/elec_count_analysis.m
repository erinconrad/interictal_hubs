function elec_count_analysis(whichPts,saved_out,out)

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
%all_surrounds = 12*[2];
main_surround = 3; %24 hour peri-revision surround
main_metric = 1;
ex_p = 1;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'rate'};
n_metrics = length(all_metrics);

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

if saved_out == 1
    
    %out = load([main_spike_results,'out.mat']);
    %out = out.out;
    
else
    out = initialize_out_struct(length(whichPts));
    
    %% Get spike details
    fprintf('Getting spike details for pt...\n');
    for i = 1:length(whichPts)
        p = whichPts(i);
        fprintf('%d of %d\n',i,length(whichPts));
        out(i) = get_gdf_details(p);
    end
    save([main_spike_results,'out'],'out');
end

%% Get names
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

%% Do main analysis

% initialize arrays
all_ov = nan(n_metrics,n_surrounds,length(whichPts),2); % 2 is pre and post rate
all_stab = nan(n_metrics,n_surrounds,length(whichPts)); % Spike stability
all_added_elecs= nan(length(whichPts),2); % 2 is all and depths
ov_stats = nan(n_metrics,n_surrounds,4); % r_all, p_all, r_depth, p_depth
stab_stats = nan(n_metrics,n_surrounds,4); % r_all, p_all, r_depth, p_depth

% Loop over metrics
for im = 1:n_metrics
    allcatlist = [];
    metric = all_metrics{im};
    
    % Loop over surround times
    for is = 1:n_surrounds

        surround = all_surrounds(is);
        
        %% Get individual patient data
        
        for i = 1:length(whichPts)
            
            switch metric
                case 'rate'
                    rate = out(i).rate;
                case 'ns'
                    
                    rate = out(i).metrics.ns;
                   
            end
            
            cblock = out(i).change_block;

            % Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            rate(ekg,:) = [];
            run_dur = out(i).run_dur;
            
            % Get surround times, starting with first non nan
            [pre,post] = get_surround_times(rate,cblock,surround);
            
            % Get the rate for all electrodes in the pre and post
            rate_pre = nanmean(rate(:,pre),2); % divide by 5 minutes to get rate per minute
            rate_post = nanmean(rate(:,post),2);
            
            % If doing rate, divide by run_dur to get spikes/minute
            if im == 1
                rate_pre = rate_pre/run_dur;
                rate_post = rate_post/run_dur;
            end
            
            % overall rate
            ov_pre = nansum(rate_pre); % sum across electrodes
            ov_post = nansum(rate_post); 
            
            % Spike stability (SRC)
            pre_post_corr = corr(rate_pre,rate_post,'Type','Spearman','rows','pairwise');

            % Fill up array with data
            all_ov(im,is,i,:) = [ov_pre ov_post];
            all_stab(im,is,i) = pre_post_corr;
            
            % Get info on added electrodes
            n_added = length(out(i).change(end).added);
            n_depths = length(out(i).change(end).added_depths);
            
            % Fill up info
            all_added_elecs(i,1) = n_added;
            all_added_elecs(i,2) = n_depths;
            
        end
        
        %% Do stats across patients
        
        % (1) Does relative spike rate increase correlate with number of electrodes
        curr_ov = squeeze(all_ov(im,is,:,:));
        ov_change = (curr_ov(:,2)-curr_ov(:,1))./(curr_ov(:,1));
        [r_ov_all,p_ov_all] = corr(ov_change,all_added_elecs(:,1));
        
        % (2) Does relative spike rate increase correlate with number of depth electrodes
        [r_ov_depth,p_ov_depth] = corr(ov_change,all_added_elecs(:,2));
        
        % (3) Does the spike stability correlate with the number of electrodes
        curr_stab = squeeze(all_stab(im,is,:));
        [r_stab_all,p_stab_all] = corr(curr_stab,all_added_elecs(:,1));
        
        % (4) Does the spike stability correlate with the number of depth electrodes
        [r_stab_depth,p_stab_depth] = corr(curr_stab,all_added_elecs(:,2));
        
        
        % Fill up arrays
        ov_stats(im,is,:) = [r_ov_all,p_ov_all,r_ov_depth,p_ov_depth];
        stab_stats(im,is,:) = [r_stab_all,p_stab_all,r_stab_depth,p_stab_depth];
        
    end
    
end

%% Initialize figure
im = main_metric;
is = main_surround;

figure
set(gcf,'position',[100 87 800 300])
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

%% Spike rate correlated with number of added electrodes
nexttile
curr_ov = squeeze(all_ov(im,is,:,:));
ov_change = (curr_ov(:,2)-curr_ov(:,1))./(curr_ov(:,1));
plot(all_added_elecs(:,1),ov_change,'o','markersize',15,'linewidth',2)
xlabel('# Added electrodes');
ylabel('Relative spike rate change')
set(gca,'fontsize',15)
xl = xlim;
yl = ylim;
text(xl(2),yl(1),sprintf('r = %1.2f, %s',...
    ov_stats(im,is,1),get_p_text(ov_stats(im,is,2))),...
    'horizontalalignment','right','fontsize',15)
pause(0.3)

%% Spike stability correlated with number of added electrodes
nexttile
curr_stab = squeeze(all_stab(im,is,:));
plot(all_added_elecs(:,1),curr_stab,'o','markersize',15,'linewidth',2)
xlabel('# Added electrodes');
ylabel('Peri-revision spike stability')
set(gca,'fontsize',15)
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('r = %1.2f, %s',...
    stab_stats(im,is,1),get_p_text(stab_stats(im,is,2))),...
    'horizontalalignment','right',...
    'verticalalignment','top','fontsize',15)

end