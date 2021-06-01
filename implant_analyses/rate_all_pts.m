function rate_all_pts(whichPts,saved_out)

%% Parameters
surround = 48;
nb = 1e4;

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

%% Get overall rate and do significance testing
for i = 1:length(whichPts)
    out(i).overall_rate = nansum(out(i).rate,1);
    out(i).nan_blocks = find(isnan(nanmean(out(i).rate,1)));
    
    out(i).overall_rate_pval = mc_overall_rate(out(i).overall_rate,...
        surround,out(i).change_block,nb);
end

%% Initialize figure
figure
set(gcf,'position',[252 547 1001 250])
tiledlayout(1,2,'TileSpacing','compact','padding','compact')

%% Example rate
nexttile
p = 1;
curr_rate = out(p).overall_rate /out(p).run_dur;
curr_times = (1:length(curr_rate)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
plot(curr_times,curr_rate,'k','linewidth',2)
hold on
nan_blocks = out(p).nan_blocks;

xlim([0 length(curr_rate)*out(p).block_dur]);
ylabel('Spikes/min')
xlabel('Hour')
title(sprintf('%Pt %d',p));
set(gca,'fontsize',20)
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(curr_rate,2))];
    bidx = bidx*out(p).block_dur;
    yl = ylim;
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),yl(2),yl(2)],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([curr_change curr_change],ylim,'r--','linewidth',4);
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')


%% Change for each patient
nexttile
all_pre = nan(length(whichPts),1);
all_post = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).overall_rate/out(i).run_dur;
    cblock = out(i).change_block;
    pre = cblock - surround:cblock -1;
    post = cblock + 1:cblock+surround;
    pre_rate = nanmean(rate(pre));
    post_rate = nanmean(rate(post));
    %{
    pr = plot([i-0.1],[pre_rate],'o','color',[0, 0.4470, 0.7410]);
    hold on
    ps = plot([i+0.1],[post_rate],'o','color',[0.8500, 0.3250, 0.098]);
    %}
    all_pre(i) = pre_rate;
    all_post(i) = post_rate;
end
%xticks([1:length(whichPts)]);
%xlabel('Patient ID')
ylabel('Spikes/min')
plot(1+0.05*rand(length(whichPts),1),all_pre,'o','color',[0, 0.4470, 0.7410],...
    'markersize',15,'linewidth',2)
hold on
plot(2+0.05*rand(length(whichPts),1),all_post,'o','color',[0.8500, 0.3250, 0.098],...
    'markersize',15,'linewidth',2)
[~,pval] = ttest(all_pre,all_post);
xlim([0.5 2.5])
yl = ylim;
ylim([yl(1) 1.10*(yl(2)-yl(1))])
plot([1 2],[yl(1) + 0.95*(yl(2)-yl(1)) yl(1) + 0.95*(yl(2)-yl(1))],'k-')
text(1.5,yl(1) + 1.02*(yl(2)-yl(1)),get_asterisks(pval,1),'horizontalalignment','center','fontsize',20)

xticks([1 2])
ylabel('Spikes/min')
xticklabels({'Pre-revision','Post-revision'})
set(gca,'fontsize',20)
%legend([pr ps],{'Pre-revision','Post-revision'},'fontsize',20);
%{
yl = ylim;
for i = 1:length(whichPts)
    pval = out(i).overall_rate_pval;
    if pval < 0.05/length(whichPts)
        text(i,yl(1)+1.05*(yl(2)-yl(1)),'*','horizontalalignment','center','fontsize',20)
    end
end
%}

%% Add subtitle labels
annotation('textbox',[0 0.93 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.5 0.93 0.1 0.1],'String','B','fontsize',30,'linestyle','none')

print(gcf,[main_spike_results,'Fig1'],'-depsc')
print(gcf,[main_spike_results,'Fig1'],'-dpng')


end