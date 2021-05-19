function pval = moran_perm_analysis(rate,change,surround,locs,spikey_idx)

nb = 1e3;
nblocks = size(rate,2);
dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));


%% Restrict to spikey
locs = locs(spikey_idx,:);
rate = rate(spikey_idx,:);

%% Calculate rate change
pre = rate(:,change-surround:change-1);
post = rate(:,change+1:change+surround);
pre_mean = nanmean(pre,2);
post_mean = nanmean(post,2);
rchange = post_mean - pre_mean;
MI = moran_index(locs,rchange,dmin);
I = MI.I;

perm_I = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    
    % recalculate
    fpre = rate(:,fchange-surround:fchange-1);
    fpost = rate(:,fchange+1:fchange+surround);
    fchange = nanmean(fpost,2) - nanmean(fpre,2);
    MI = moran_index(locs,fchange,dmin);
    perm_I(ib) = MI.I;
    
end

perm_I = sort(perm_I);

if 0
    plot(perm_I,'o')
    hold on
    plot(xlim,[I I])
    
end

num_as_sig = sum(abs(perm_I) >= abs(I)); % those more positive or more negative
pval = (num_as_sig+1)/(nb+1);

end