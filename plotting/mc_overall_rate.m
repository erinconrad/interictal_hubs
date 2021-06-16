function pval = mc_overall_rate(rate,surround,change,nb,only_pre)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
[pre,post] = get_surround_times(rate,change,surround);
pre = rate(pre);
post = rate(post);

pre_mean = nanmean(pre);
post_mean = nanmean(post);
true_diff = post_mean - pre_mean;

perm_diff = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    if only_pre
        fchange = randi([surround+1,change-1]);
    else
        fchange = randi([surround+1,nblocks-surround]);
    end
    [fpre,fpost] = get_surround_times(rate,fchange,surround);
    
    fpre = rate(fpre);
    fpost = rate(fpost);
    
    
    perm_diff(ib) = nanmean(fpost) - nanmean(fpre);
    
end

%% sort the permutation results
[sorted_perm_diff] = sort(perm_diff);

% How many differences as big as true diff (negative or positive, two
% tailed)
num_as_sig = sum(abs(sorted_perm_diff) >= abs(true_diff));
pval = (num_as_sig+1)/(nb+1);


end