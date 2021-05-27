function show_overall_rate(all_rate,block_dur,change_block,run_dur,name,results_folder,surround)


%% Significance testing
pval = rate_analysis(all_rate,change_block,surround);
%pval = overall_pre_post_rate(all_rate,change_block);

%{
pre = 1:change_block-1;
post = change_block+1:size(all_rate,2);
rate_pre = nanmean(all_rate(:,pre),1);
rate_post = nanmean(all_rate(:,post),1);

pval = ranksum(rate_pre,rate_post);
%}

figure
set(gcf,'position',[172 233 1181 423])
times = 1:size(all_rate,2);
times = times*block_dur;
%plot(times,all_rate/run_dur,'--','color',[0, 0.4470, 0.7410])
hold on
plot(times,nansum(all_rate,1)/run_dur,'k','linewidth',2)

hold on

nan_blocks = find(isnan(nanmean(all_rate,1)));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(all_rate,2))];
    bidx = bidx*block_dur;
    yl = ylim;
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),yl(2),yl(2)],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',4);

xlim([0 size(all_rate,2)*block_dur]);
ylabel('Spikes/min')
xlabel('Hour')
title(sprintf('%s p = %1.3f',name,pval));
set(gca,'fontsize',20)
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')


output_folder = [results_folder,'overall_rate/'];
if exist(output_folder,'dir') == 0
    mkdir(output_folder)
end

print(gcf,[output_folder,name,'_rate_',sprintf('%d',surround)],'-dpng')
close(gcf)

end

function pval = rate_analysis(rate,change,surround)



do_vecnorm = 0;

nb = 1e3;

nblocks = size(rate,2);

%% Compare spike rate pre and post change
pre = rate(:,change-surround:change-1);
post = rate(:,change+1:change+surround);

if do_vecnorm
    a = nanmean(post,2) - nanmean(pre,2);
    true_diff = vecnorm(a(~isnan(a)));
else
    pre_mean = nanmean(pre(:));
    post_mean = nanmean(post(:));
    true_diff = post_mean - pre_mean;
end

perm_diff = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    
    % recalculate
    fpre = rate(:,fchange-surround:fchange-1);
    fpost = rate(:,fchange+1:fchange+surround);
    
    if do_vecnorm
        a = nanmean(fpost,2) - nanmean(fpre,2);
        perm_diff(ib) = vecnorm(a(~isnan(a)));
    else
        perm_diff(ib) = nanmean(fpost(:)) - nanmean(fpre(:));
    end
    
end

%% sort the permutation results
[sorted_perm_diff] = sort(perm_diff);

if 0
    plot(sorted_perm_diff,'o')
    hold on
    plot(xlim,[true_diff true_diff])
    
end

if do_vecnorm
    % one tailed test, is the diff larger (diff > 0)
    num_as_sig = sum(sorted_perm_diff>=true_diff);
    
else
    % For a two tailed test, find the number above true_diff or <-true_diff
    num_as_sig = sum(abs(sorted_perm_diff) >= abs(true_diff));
    %{
    if true_diff >= 0
        num_as_sig = sum(sorted_perm_diff>=true_diff | sorted_perm_diff <= -true_diff);
    else
        num_as_sig = sum(sorted_perm_diff<=true_diff | sorted_perm_diff >= -true_diff);
    end
    %}
end

pval = (num_as_sig+1)/(nb+1);

end