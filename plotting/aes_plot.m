function aes_plot(rate,block_dur,change,run_dur,unchanged_locs,added_locs,...
    name,results_folder,labels)

%% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(labels);
rate(ekg,:) = [];
unchanged_locs(ekg,:) = [];
added_locs(ekg,:) = [];

figure
set(gcf,'position',[100 100 1300 600]) 
t = tiledlayout(2,2,'TileSpacing','compact','Padding','Compact');

%% Elec position
nexttile(1)
pu = scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),...
    150,'linewidth',2);
hold on
pa = scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),150,'p',...
    'linewidth',2);
set(gca,'Visible','off')
title('Electrode locations')
legend('Pre-existing','Added','fontsize',20,'location','southwest')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gca,'fontsize',20)

%% Rate
nexttile(2)
times = 1:size(rate,2);
times = times*block_dur;
plot(times,nansum(rate,1)/run_dur,'k','linewidth',2)

hold on
nan_blocks = find(isnan(nanmean(rate,1)));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(rate,2))];
    bidx = bidx*block_dur;
    yl = ylim;
    %{
    ap = area(bidx,ylim,'basevalue',0,...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
   %}
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),yl(2),yl(2)],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
    
end
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',4);

xlim([0 size(rate,2)*block_dur]);
ylabel('Spikes/min')
xlabel('Hour')
set(gca,'fontsize',20)
title('Overall spike rate')
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')
ylim([0 max(nansum(rate,1)/run_dur)])

%% Raster
nexttile(3,[1 2])
times = 1:size(rate,2);
times = times*block_dur;
h = turn_nans_white(rate/run_dur);
set(h,'XData',[0:size(rate,2)*block_dur])
xlim([0 size(rate,2)*block_dur])
hold on
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',5);
yticklabels([])
xlabel('Hour')
ylabel('Electrode')
title('Spike rate by electrode')
c = colorbar;
ylabel(c,'Spikes/min','fontsize',20)

set(gca,'fontsize',20)

%% Add labels
annotation('textbox',[0.04 0.9 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.48 0.9 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
annotation('textbox',[0.04 0.4 0.1 0.1],'String','C','fontsize',30,'linestyle','none')


outfolder = [results_folder,'aes/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder);
end
print(gcf,[outfolder,name,'_aes'],'-depsc')


end