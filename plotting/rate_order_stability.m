function rate_order_stability(all_rate,all_rl,change_block,surround,block_dur,name,results_folder)

% Significance testing
pval = compare_rhos(all_rate,change_block,surround);

% times
times = 1:size(all_rate,2);
times = times * block_dur;

% Take the mean across the pre-change surround
mean_rate = nanmean(all_rate,2);
mean_rl = nanmean(all_rl,2);

% For each block, take the correlation between the mean rate and
% the current rate
rho = corr(all_rate,mean_rate,'rows','pairwise');
rho_rl = corr(all_rl,mean_rl,'rows','pairwise');

figure
set(gcf,'position',[202 443 1239 300]);
%subplot(2,1,1)
plot(times,rho)
hold on
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);
ylabel('Rate order stability')
xlabel('Hour')
title(sprintf('Mean ros: %1.2f, p = %1.3f',nanmean(rho),pval))
set(gca,'fontsize',20)
%{
subplot(2,1,2)
plot(times,rho_rl)
hold on
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);
ylabel('Sequence stability')
xlabel('Hour')
title(sprintf('Mean ss: %1.2f',nanmean(rho_rl)))
set(gca,'fontsize',20)
%}

out_folder = [results_folder,'rate_order/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

print(gcf,[out_folder,name,'_rate_order'],'-dpng');
close(gcf);


end


