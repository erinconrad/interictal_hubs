function pval = mc_overall_rate(rate,surround,change,nb,full_rate)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
[pre,post] = get_surround_times(full_rate,change,surround);
pre = rate(pre);
post = rate(post);

pre_mean = nanmean(pre);
post_mean = nanmean(post);
true_diff = post_mean - pre_mean;

perm_diff = nan(nb,1);

for ib = 1:nb
    
    while 1
        % Make a fake change time
        fchange = randi([surround+1,nblocks-surround]);

        [fpre,fpost] = get_surround_times(full_rate,fchange,surround);

        if isnan(fpre), continue; end
        
        fpre = rate(fpre);
        fpost = rate(fpost);

        perm_diff_temp = nanmean(fpost) - nanmean(fpre);
        
        if isnan(perm_diff_temp)
            continue
        else
            break
        end
    end
    perm_diff(ib) = perm_diff_temp;
    
end

%% sort the permutation results
[sorted_perm_diff] = sort(perm_diff);

if sum(isnan(perm_diff))>0
    error('why');
end

% How many differences as big as true diff (negative or positive, two
% tailed)
num_as_sig = sum(abs(sorted_perm_diff) >= abs(true_diff));
pval = (num_as_sig+1)/(nb+1);


end