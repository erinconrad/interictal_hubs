function ros_all_pts(whichPts,saved_out)

%{
Generates Figure 2 for implant effect paper -> raster plot for a single
patient and rate order stability for all patients and stability of spikiest
electrode for all patients
%}

%% Parameters
surround = 48;
nb = 1e4;
pl_surround = 96;

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
    
    out = load([main_spike_results,'out.mat']);
    out = out.out;
    
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

%% Initialize figure
figure
set(gcf,'position',[100 100 1001 550])
tiledlayout(2,2,'TileSpacing','compact','padding','compact')

%% Example spike rate raster
nexttile
p = 1;
rate_raster = out(p).rate./out(p).run_dur;
% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(out(p).unchanged_labels);
rate_raster(ekg,:) = [];
curr_times = (1:size(rate_raster,2)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
h = turn_nans_white(rate_raster);
set(h,'XData',[0:curr_times(end)]);
xlim([0 curr_times(end)])
hold on
cp = plot([curr_change curr_change],ylim,'r--','linewidth',3);
xlabel('Hour')
legend(cp,'Revision','fontsize',20,'location','southeast')
yticklabels([])
ylabel('Electrode')
set(gca,'fontsize',20)
c = colorbar;
ylabel(c,'Spikes/min','fontsize',20)

%% Example node strength raster
nexttile

%% Rate order stability
nexttile
all_ros = nan(length(whichPts),pl_surround*2+1);
all_ps = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    
    % Take mean over pre-revision period
    first_rate = nanmean(rate(:,cblock-pl_surround:cblock-1),2);
    nblocks = pl_surround*2+1;
       
    % Get rate order stability
    ros = nan(nblocks,1);
    hcount = 0;
    for h = cblock-pl_surround:cblock+pl_surround
        hcount = hcount + 1;
        ros(hcount) = corr(first_rate,rate(:,h),'Type','Spearman','rows','pairwise');
    end
    
    % put it in the cell
    all_ros(i,:) = ros;
    
    % do a ros permutation test to see whether the pre-post change is
    % larger than expected for randomly chosen times
    pval_curr = compare_rhos(rate,cblock,surround,nb);
    all_ps(i) = pval_curr;
    
end

% Fisher test to combine pvalues
X_2 = -2 * sum(log(all_ps));
sum_p = 1-chi2cdf(X_2,2*length(all_ps));

% Plot the times
times = -pl_surround:pl_surround;
m = nanmean(all_ros,1);
st = nanstd(all_ros,[],1);
[mp,stp] = shaded_error_bars(times,m,st,[]);
hold on
set(gca,'fontsize',20)
xlabel('Hours surrounding revision')
ylabel('Rate order stability')
yl = get(gca,'ylim');
ylim([yl(1) yl(1) + 1.11*(yl(2)-yl(1))])
plot([0 0],[yl(1) 0.98*(yl(2)-yl(1))],'r--','linewidth',3)
plot([-surround surround],[1.01*(yl(2)-yl(1)) 1.01*(yl(2)-yl(1))],'k-','linewidth',2)
text(0,1.03*(yl(2)-yl(1)),get_asterisks(sum_p,1),...
    'horizontalalignment','center','fontsize',30)

xlim([times(1) times(end)])

%% NS order stability
nexttile

%% Spikiest electrode stability
%{
nexttile
nblocks = pl_surround*2+1;
agree_all = nan(length(whichPts),nblocks);
pval_all = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    
    % For each time period, get the spikiest electrode
    [~,spikiest] = max(rate,[],1);
    
    % Get mean rate in pre-revision period
    
    mean_pre_rate = nanmean(rate(:,cblock-pl_surround:cblock-1),2);
    [~,spikiest_pre] = max(mean_pre_rate,[],1);
    
    % For each time period, get whether spikiest electrode agrees with pre
    agree = nan(nblocks,1);
    hcount = 0;
    for h = cblock-pl_surround:cblock+pl_surround
        hcount = hcount + 1;
        agree(hcount) = isequal(spikiest_pre,spikiest(h));
    end
    
    agree_all(i,:) = agree;
    
    % Significance testing
    pval = compare_agree(rate,cblock,surround,nb);
    pval_all(i) = pval;
    
end

% Fisher
X_2 = -2 * sum(log(pval_all));
sum_p = 1-chi2cdf(X_2,2*length(pval_all));

% Plot the agree
times = -pl_surround:pl_surround;
plot(times,nanmean(agree_all,1),'linewidth',3)
hold on
set(gca,'fontsize',20)
xlabel('Hours surrounding revision')
ylabel('Spike hub stability')
yl = get(gca,'ylim');
ylim([yl(1) yl(1) + 1.11*(yl(2)-yl(1))])
plot([0 0],[yl(1) 0.98*(yl(2)-yl(1))],'r--','linewidth',3)
plot([-surround surround],[1.01*(yl(2)-yl(1)) 1.01*(yl(2)-yl(1))],'k-','linewidth',2)
text(0,1.03*(yl(2)-yl(1)),get_asterisks(sum_p,1),...
    'horizontalalignment','center','fontsize',30)

xlim([times(1) times(end)])

%}

end