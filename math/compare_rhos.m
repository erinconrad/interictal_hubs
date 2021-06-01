function pval = compare_rhos(rate,change,surround,nb)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
pre = rate(:,change-surround:change-1);
post = rate(:,change+1:change+surround);
pre_mean = nanmean(pre,2);
post_mean = nanmean(post,2);
true_rho = corr(pre_mean,post_mean,'rows','pairwise');

perm_rho = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    %fchange = randi([change,nblocks-surround]);
    
    % recalculate
    fpre = rate(:,fchange-surround:fchange-1);
    fpost = rate(:,fchange+1:fchange+surround);
    
    perm_rho(ib) = corr(nanmean(fpre,2),nanmean(fpost,2),'rows','pairwise');
    
end

%% sort the permutation results
[sorted_perm_rho] = sort(perm_rho);

if 0
    plot(sorted_perm_rho,'o')
    hold on
    plot(xlim,[true_rho true_rho])
end

% this is just a one-tailed test which I think is fair
num_as_sig = sum(sorted_perm_rho<=true_rho);
pval = (num_as_sig+1)/(nb+1);

end