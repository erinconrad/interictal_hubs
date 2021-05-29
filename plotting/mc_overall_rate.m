function pval = mc_overall_rate(rate,surround,change,nb)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
pre = rate(change-surround:change-1);
post = rate(change+1:change+surround);

pre_mean = nanmean(pre);
post_mean = nanmean(post);
true_diff = post_mean - pre_mean;

perm_diff = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    
    % recalculate
    fpre = rate(fchange-surround:fchange-1);
    fpost = rate(fchange+1:fchange+surround);
    
    perm_diff(ib) = nanmean(fpost) - nanmean(fpre);
    
end

%% sort the permutation results
[sorted_perm_diff] = sort(perm_diff);
num_as_sig = sum(abs(sorted_perm_diff) >= abs(true_diff));
pval = (num_as_sig+1)/(nb+1);


end