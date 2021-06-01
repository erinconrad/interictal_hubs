function cosi_over_time(rate,cosi,spikey_idx,labels,name,results_folder,...
    block_dur,change_block,run_dur,surround)

thresh = 0.8;

nan_blocks = find(isnan(nanmean(rate,1)));

%% Just do spikey
rate = rate(spikey_idx,:);
cosi = cosi(spikey_idx);
labels = labels(spikey_idx);


%% Further restrict to just those with a high cospike index
hc = cosi > thresh;

if sum(hc) == 0
    fprintf('\nSkipping %s as no cospikey channels\n',name);
    return
end

rate = rate(hc,:);
labels = labels(hc);

%% Plot
figure
set(gcf,'position',[172 233 1181 423])
times = 1:size(rate,2);
times = times*block_dur;
plot(times,rate/run_dur,'--','color',[0, 0.4470, 0.7410])
hold on
plot(times,nanmean(rate,1)/run_dur,'k','linewidth',2)
%{
text(size(rate,2)*block_dur,rate(:,end)/run_dur,labels,'fontsize',20)

%}
xlim([0 size(rate,2)*block_dur]+10);

for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(rate,2))];
    bidx = bidx*block_dur;
    ap = area(bidx,ylim,'basevalue',0,...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);


ylabel('Spikes/min')
xlabel('Hour')
set(gca,'fontsize',20)

out_folder = [results_folder,'time_cosi/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end
print(gcf,[out_folder,name],'-dpng')
close(gcf)


end