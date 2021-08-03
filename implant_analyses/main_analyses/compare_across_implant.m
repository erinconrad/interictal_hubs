function compare_across_implant(whichPts,saved_out,out)

%% User change parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
%all_surrounds = 12*[2];
main_surround = 3; %24 hour peri-revision surround
main_metric = 1;
ex_p = 1;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'rate','ns'};
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

%% Do analyses

% Initialize output
n_patients = length(whichPts);
all_rates = nan(n_metrics,n_surrounds,n_patients,4); % 4 is start, pre, post, end
all_rhos = nan(n_metrics,n_surrounds,n_patients,2); %2 is start-pre, post-end
within_implant_rate_stats = nan(n_metrics,n_surrounds,2,3); %2 is first vs second implant, 3 is p-value, tstat, df
between_implant_rate_stats = nan(n_metrics,n_surrounds,3); %3 is p-value, tstat, df
all_rhos_stats = nan(n_metrics,n_surrounds,3); %3 is p-value, tstat, df


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
            
            
            % Get pre, post, start, ending times
            [start,pre,post,ending] = across_implant_surround(rate,cblock,surround);
            
            % Get mean rates (across time periods) in these times
            rate_start = nanmean(rate(:,start),2)/run_dur;
            rate_pre = nanmean(rate(:,pre),2)/run_dur;
            rate_post = nanmean(rate(:,post),2)/run_dur;
            rate_end = nanmean(rate(:,ending),2)/run_dur;
            % This is the average spike rate in each channel for a 5 minute
            % run
            
            % Get the mean overall rate in each
            overall_rate_start = nansum(rate_start);
            overall_rate_pre = nansum(rate_pre);
            overall_rate_post = nansum(rate_post);
            overall_rate_end = nansum(rate_end);
            % And so, this would be in units of spikes/minute
            
            % Get the SRC between start and pre, and post and end
            start_pre_corr = corr(rate_start,rate_pre,'Type','Spearman','rows','pairwise');
            post_end_corr = corr(rate_post,rate_end,'Type','Spearman','rows','pairwise');
            
            % Fill up matrix
            all_rates(im,is,i,:) = [overall_rate_start,overall_rate_pre,...
                overall_rate_post,overall_rate_end];
            all_rhos(im,is,i,:) = [start_pre_corr,post_end_corr];
            
        end
        
        %% Do stats across patients
        
        % Compare rates within first implant
        all_first_rates = squeeze(all_rates(im,is,:,1:2));
        [~,p,~,stats] = ttest(all_first_rates(:,1),all_first_rates(:,2));
        
        % fill matrix
        within_implant_rate_stats(im,is,1,:) = [p,stats.tstat,stats.df];
        
        % compare rates within 2nd implant
        all_second_rates = squeeze(all_rates(im,is,:,3:4));
        [~,p,~,stats] = ttest(all_second_rates(:,1),all_second_rates(:,2));
        
        % fill matrix
        within_implant_rate_stats(im,is,2,:) = [p,stats.tstat,stats.df];
        
        % Compare rate from start to end
        [~,p,~,stats] = ttest(all_first_rates(:,1),all_second_rates(:,2));
        
        % fill matrix
        between_implant_rate_stats(im,is,:) = [p,stats.tstat,stats.df];
        
        % Compare rhos between first and second implant
        curr_all_rhos = squeeze(all_rhos(im,is,:,:));
        curr_all_zs = fisher_transform(curr_all_rhos,nan); % convert rhos to z's
        [~,p,~,stats] = ttest(curr_all_zs(:,1),curr_all_zs(:,2));
        all_rhos_stats(im,is,:) = [p,stats.tstat,stats.df];
        
    end
    
end

%% Specify metric and surround for figures
im = main_metric;
is = main_surround;

%% initialize figure
figure
set(gcf,'position',[100 87 733 900])
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

% colors
cols = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980];

%% example raster
nexttile
rate = out(ex_p).rate./out(ex_p).run_dur;
ekg = identify_ekg_scalp(out(ex_p).unchanged_labels);
rate(ekg,:) = [];
curr_times = (1:size(rate,2)) * out(ex_p).block_dur;
curr_change = out(ex_p).change_block*out(ex_p).block_dur;
h = turn_nans_white(rate);
set(h,'XData',[0:curr_times(end)]);
xlim([0 curr_times(end)])
hold on
yl = get(gca,'ylim');
cp = plot([curr_change curr_change],yl,'r--','linewidth',3);
xlabel('Hour')
yticklabels([])
if im == 1
    ylabel('Electrode')
end
set(gca,'fontsize',15)
c = colorbar;
ylabel(c,'Spikes/min','fontsize',15);


%% Rate analysis
nexttile

% plot the data
curr_rates = squeeze(all_rates(im,is,:,:));
pfirst = plot(([1 2].*ones(n_patients,1))',([curr_rates(:,1) curr_rates(:,2)])',...
    'color',cols(1,:),'linewidth',2);
hold on
psecond = plot(([3 4].*ones(n_patients,1))',([curr_rates(:,3) curr_rates(:,4)])',...
    'color',cols(2,:),'linewidth',2);
xlim([0 5])

% Get ylim stuff
yl = get(gca,'ylim');
new_yl = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
bar_y = yl(1) + 1.1*(yl(2)-yl(1));
p_y = yl(1) + 1.15*(yl(2)-yl(1));
higher_bar = yl(1) + 1.2*(yl(2)-yl(1));
higher_p = yl(1) + 1.25*(yl(2)-yl(1));
set(gca,'ylim',new_yl);

% plot the stats
plot([1 2],[bar_y bar_y],'k','linewidth',1)
text(1.5,p_y,get_asterisks(squeeze(within_implant_rate_stats(im,is,1,1)),1),...
    'horizontalalignment','center','fontsize',15);

plot([3 4],[bar_y bar_y],'k','linewidth',1)
text(3.5,p_y,get_asterisks(squeeze(within_implant_rate_stats(im,is,2,1)),1),...
    'horizontalalignment','center','fontsize',15);

plot([1 4],[higher_bar higher_bar],'k','linewidth',1)
text(2.5,higher_p,get_asterisks(squeeze(between_implant_rate_stats(im,is,1)),1),...
    'horizontalalignment','center','fontsize',15);

% Labels
xticks([1 2 3 4])
xticklabels({'Early','Late','Early','Late'})
legend([pfirst(1),psecond(1)],{'Implant 1','Implant 2'},'fontsize',15)
ylabel('Spikes/min')
set(gca,'fontsize',15);


%% ROS analysis
nexttile
curr_rhos = squeeze(all_rhos(im,is,:,:));
plot(1+0.05*rand(n_patients,1),curr_rhos(:,1),'o','markersize',15,'linewidth',2)
hold on
plot(2+0.05*rand(n_patients,1),curr_rhos(:,2),'o','markersize',15,'linewidth',2)
xlim([0 3])

% Get ylim stuff
set(gca,'ylim',[0 1]);
yl = ylim;
new_yl = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
bar_y = yl(1) + 1.1*(yl(2)-yl(1));
p_y = yl(1) + 1.15*(yl(2)-yl(1));
ylim(new_yl)

% Plot the stats
plot([1 2],[bar_y bar_y],'k','linewidth',1)
text(1.5,p_y,get_asterisks(squeeze(all_rhos_stats(im,is,1)),1),...
    'horizontalalignment','center','fontsize',15);

% Labels
ylabel('Early-to-late spike stability')
xticks([1 2])
xticklabels({'Implant 1','Implant 2'})
set(gca,'fontsize',15);

%% Make supplemental tables
im1_spike_T = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),...
    within_implant_rate_stats(im,:,1,3)',within_implant_rate_stats(im,:,1,2)',...
    within_implant_rate_stats(im,:,1,1)','UniformOutput',false));
im2_spike_T = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),...
    within_implant_rate_stats(im,:,2,3)',within_implant_rate_stats(im,:,2,2)',...
    within_implant_rate_stats(im,:,2,1)','UniformOutput',false));
bet_rate_spike_T = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),...
    between_implant_rate_stats(im,:,3)',between_implant_rate_stats(im,:,2)',...
    between_implant_rate_stats(im,:,1)','UniformOutput',false));
rho_spike_T = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),...
    all_rhos_stats(im,:,3)',all_rhos_stats(im,:,2)',...
    all_rhos_stats(im,:,1)','UniformOutput',false));

writetable(im1_spike_T,[main_spike_results,'Supplemental Table 2.xlsx'],'Range','B2:B12','WriteVariableNames',false)
writetable(im2_spike_T,[main_spike_results,'Supplemental Table 2.xlsx'],'Range','C2:C12','WriteVariableNames',false)
writetable(bet_rate_spike_T,[main_spike_results,'Supplemental Table 2.xlsx'],'Range','D2:D12','WriteVariableNames',false)
writetable(rho_spike_T,[main_spike_results,'Supplemental Table 2.xlsx'],'Range','E2:E12','WriteVariableNames',false)


end