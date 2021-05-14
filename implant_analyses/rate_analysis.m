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
    if true_diff >= 0
        num_as_sig = sum(sorted_perm_diff>=true_diff | sorted_perm_diff <= -true_diff);
    else
        num_as_sig = sum(sorted_perm_diff<=true_diff | sorted_perm_diff >= -true_diff);
    end
end

pval = (num_as_sig+1)/(nb+1);

end