function change_dist_time(all_dist,block_dur,change_block,name,results_folder,surround,all_rate)

do_mean = 1;
if do_mean
    mean_dist = cellfun(@(x) nanmean(x),all_dist);
    mtext = '_mean';
else
    mean_dist = cellfun(@(x) mode(x),all_dist);
    mtext = '_mode';
end

pval = stats_time_comp(mean_dist,change_block,surround);

figure
set(gcf,'position',[172 233 1181 423])
times = 1:size(all_rate,2);
times = times*block_dur;
plot(times,mean_dist)
hold on

nan_blocks = find(isnan(nanmean(all_rate,1)));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(all_rate,2))];
    bidx = bidx*block_dur;
    ap = area(bidx,ylim,'basevalue',0,...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);

xlim([0 size(all_rate,2)*block_dur]);
xlabel('Hour')
ylabel('Mean distance from nearest added elecs')
title(sprintf('%s p = %1.3f',name,pval));

set(gca,'fontsize',20)
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')

output_folder = [results_folder,'dist_time/'];
if exist(output_folder,'dir') == 0
    mkdir(output_folder)
end

print(gcf,[output_folder,name],'-dpng')
close(gcf)

end