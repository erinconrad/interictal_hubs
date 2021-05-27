function pval = compare_agree(rate,change,surround)

nb = 1e3;
nblocks = size(rate,2);

[~,spikiest] = max(rate,[],1);

%% turn spikiest for nans periods into nans
nantimes = isnan(nanmean(rate,1));
spikiest(nantimes) = nan;

%% Get the spikiest electrode pre-revision
pre = change-surround:change-1;
post = change+1:change+surround;
[~,spikiest_pre] = max(nanmean(rate(:,pre),2));

%% For each time period, get the spikiest electrode and whether it agrees with the pre-revision
agree = spikiest == spikiest_pre;
true_agree_pre = nanmean(agree(pre));
true_agree_post = nanmean(agree(post));

true_agree_diff = true_agree_pre-true_agree_post;

perm_agree_diff = nan(nb,1);
for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    %fchange = randi([change,nblocks-surround]);
    
    % recalculate
    fpre = fchange-surround:fchange-1;
    fpost = fchange+1:fchange+surround;
    
    [~,fspikiest_pre] = max(nanmean(rate(:,fpre),2));
    agree = spikiest == fspikiest_pre;
    agree_pre = nanmean(agree(fpre));
    agree_post = nanmean(agree(fpost));
    
    perm_agree_diff(ib) = agree_pre-agree_post;
    
end

%% sort the permutation results
perm_agree_diff = sort(perm_agree_diff);

%% this is just a one-tailed test (only really meaningfully expect pre agreement to be higher than post)
num_as_sig = sum(perm_agree_diff>=true_agree_diff);
pval = (num_as_sig+1)/(nb+1);


if 0
    plot(perm_agree_diff,'o')
    hold on
    plot(xlim,[true_agree_diff true_agree_diff])
    title(sprintf('%1.3f',pval))
end

end
