function rate_all_pts(whichPts,saved_out)

%{
This makes the rate figure and analysis
%}

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10]; % number of half hour segments surrounding implant
%all_surrounds = 12*2;
main_surround = 3;
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
surround = all_surrounds(main_surround);
for i = 1:length(whichPts)
    out(i).overall_rate = nansum(out(i).rate,1); % sum across electrodes to get total number of spikes
    out(i).nan_blocks = find(isnan(nanmean(out(i).rate,1)));
    
    out(i).overall_rate_pval = mc_overall_rate(out(i).overall_rate,...
        surround,out(i).change_block,nb);
end

%% Initialize figure
figure
set(gcf,'position',[50 547 1200 600])
tiledlayout(4,3,'TileSpacing','compact','padding','compact')
tile_order = [1 2 4 5 7 8 10 11 3 6 9];

%% Example rates
for p = 1:length(out)
nexttile(tile_order(p))
curr_rate = out(p).overall_rate /out(p).run_dur;
curr_times = (1:length(curr_rate)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
plot(curr_times,curr_rate,'k','linewidth',1)
hold on
nan_blocks = out(p).nan_blocks;
[pre,post] = get_surround_times(out(p).overall_rate,out(p).change_block,surround);
pre = pre*out(p).block_dur;
post = post*out(p).block_dur;
xlim([0 length(curr_rate)*out(p).block_dur]);
%if ismember(tile_order(p) ,[1 4 7 10])
    ylabel('Spikes/min')
%end
%if ismember(tile_order(p),[10 11])
    xlabel('Hour')
%end

set(gca,'fontsize',15)
yl = ylim;
new_yl = [yl(1) 1.15*(yl(2)-yl(1))];
ylim(new_yl);
yl = ylim;
top = yl(1) + 0.75*(yl(2)-yl(1));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(curr_rate,2))];
    bidx = bidx*out(p).block_dur;
    
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),top,top],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([curr_change curr_change],[yl(1) top],'r--','linewidth',4);
plot([pre(1) post(end)],[yl(1) + 0.8*(yl(2)-yl(1)) yl(1) + 0.8*(yl(2)-yl(1))],...
    'k-','linewidth',2);
ast = get_asterisks(out(p).overall_rate_pval,length(out));
text(curr_change,yl(1) + 0.9*(yl(2)-yl(1)),ast,...
    'horizontalalignment','center','fontsize',15)
%legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('Patient %d',p),'fontsize',15,'VerticalAlignment','Top')
end

%% Rate change for all patients
nexttile(tile_order(end),[2 1])
all_pre = nan(length(whichPts),1);
all_post = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).overall_rate/out(i).run_dur;
    cblock = out(i).change_block;

    % Get surround times, starting with first non nan
    [pre,post] = get_surround_times(rate,cblock,surround);
    
    pre_rate = nanmean(rate(pre));
    post_rate = nanmean(rate(post));

    all_pre(i) = pre_rate;
    all_post(i) = post_rate;
end

ylabel('Spikes/min')
plot(1+0.05*rand(length(whichPts),1),all_pre,'o','color',[0, 0.4470, 0.7410],...
    'markersize',15,'linewidth',2)
hold on
plot(2+0.05*rand(length(whichPts),1),all_post,'o','color',[0.8500, 0.3250, 0.098],...
    'markersize',15,'linewidth',2)
[~,pval] = ttest(all_pre,all_post);
xlim([0.5 2.5])
yl = ylim;
ylim([yl(1) 1.15*(yl(2)-yl(1))])
plot([1 2],[yl(1) + 1*(yl(2)-yl(1)) yl(1) + 1*(yl(2)-yl(1))],'k-','linewidth',2)
text(1.5,yl(1) + 1.07*(yl(2)-yl(1)),get_asterisks(pval,1),'horizontalalignment','center','fontsize',15)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),'All patients','verticalalignment','top','fontsize',15)

xticks([1 2])
ylabel('Spikes/min')
xticklabels({'Pre-revision','Post-revision'})
set(gca,'fontsize',15)


%% Add subtitle labels
%annotation('textbox',[0 0.93 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
%annotation('textbox',[0.5 0.93 0.1 0.1],'String','B','fontsize',30,'linestyle','none')

print(gcf,[main_spike_results,'Fig1'],'-depsc')
print(gcf,[main_spike_results,'Fig1'],'-dpng')

%% Do a t-test for each surround so that I can get a table of p-values
all_ps = nan(length(all_surrounds),7);
for s = 1:length(all_surrounds)
    surround = all_surrounds(s);
    all_pre = nan(length(whichPts),1);
    all_post = nan(length(whichPts),1);
    for i = 1:length(whichPts)
        rate = out(i).overall_rate/out(i).run_dur;
        cblock = out(i).change_block;
        
        % Get surround times, starting with first non nan
        [pre,post] = get_surround_times(rate,cblock,surround);

        pre_rate = nanmean(rate(pre));
        post_rate = nanmean(rate(post));
        
        all_pre(i) = pre_rate;
        all_post(i) = post_rate;
    end
    [~,pval,~,stats] = ttest(all_pre,all_post);
    tstat = stats.tstat;
    df = stats.df;
    all_ps(s,:) = [pval tstat df mean(all_pre) std(all_pre) mean(all_post) std(all_post)];
end

%% Save table of p-values and tstats
ptext = arrayfun(@(x) sprintf('%1.3f',x),all_ps(:,1),...
    'UniformOutput',false);
ttext = arrayfun(@(x) sprintf('%1.2f',x),all_ps(:,2),...
    'UniformOutput',false);
all_text = arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),all_ps(:,3),all_ps(:,2), all_ps(:,1),...
    'UniformOutput',false);
T = table(ttext,ptext,...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'VariableNames',{'t','p'});
all_T = table(all_text,...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));

writetable(T,[main_spike_results,'rate.csv'],'WriteRowNames',true)  

%% Also save to supplemental Table xls
writetable(all_T,[main_spike_results,'Supplemental Table 1.xlsx'],'Range','B2:B12','WriteVariableNames',false)

%% Results Sentence
fprintf(['Aggregated across patients, there was no difference in the'...
    ' pre- (M = %1.1f, SD = %1.1f spikes/min) and post-revision'...
    ' (M = %1.1f, SD = %1.1f spikes/min) spike rate (paired t-test,'...
    ' t(%d) = %1.1f, p = %1.2f)\n'],...
    all_ps(1,4),all_ps(1,5),all_ps(1,6),all_ps(1,7),...
    all_ps(1,3),all_ps(1,2),all_ps(1,1))
end