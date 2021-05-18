function pval = stats_time_comp(thing,change,surround)

nb = 1e3;
nblocks = length(thing);

%% Compare pre and post
pre = nanmean(thing(change-surround:change-1));
post = nanmean(thing(change+1:change+surround));
true_diff = post - pre;

perm_diff = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    
    % recalculate
    fpre = nanmean(thing(fchange-surround:fchange-1));
    fpost = nanmean(thing(fchange+1:fchange+surround));
    perm_diff(ib) = nanmean(fpost(:)) - nanmean(fpre(:));
    
    
end

%% sort the permutation results
[sorted_perm_diff] = sort(perm_diff);


if 0
    figure
    plot(sorted_perm_diff,'o')
    hold on
    plot(xlim,[true_diff true_diff])
    
end

num_as_sig = sum(sorted_perm_diff>=abs(true_diff) | sorted_perm_diff <= -abs(true_diff));
pval = (num_as_sig+1)/(nb+1);

end