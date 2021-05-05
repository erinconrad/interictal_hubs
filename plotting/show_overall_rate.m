function show_overall_rate(all_rate,block_dur,last_block,change_block,run_dur,name,results_folder)

figure
set(gcf,'position',[172 233 1181 423])
times = 1:size(all_rate,2);
times = times*block_dur;
plot(times,nanmean(all_rate,1)/run_dur,'linewidth',1)

hold on
%{
for b = 1:length(last_block)
    fp = plot([last_block(b)*block_dur last_block(b)*block_dur],ylim,'k--','linewidth',3);
end
%}

nan_blocks = find(isnan(nanmean(all_rate,1)));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(all_rate,2))];
    bidx = bidx*block_dur;
    ap = area(bidx,ylim,'basevalue',0,...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);

xlim([0 size(all_rate,2)*block_dur]);
ylabel('Spikes/min')
xlabel('Hour')
title(sprintf('%s',name));
set(gca,'fontsize',20)
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')

output_folder = [results_folder,'overall_rate/'];
if exist(output_folder,'dir') == 0
    mkdir(output_folder)
end

print(gcf,[output_folder,name,'_rate'],'-dpng')

end