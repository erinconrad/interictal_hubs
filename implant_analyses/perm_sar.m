function perm_sar(thing,rate,change,surround,locs,spikey_idx,outfolder,...
    name,do_rel,ttext)

dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));

thing = thing(spikey_idx);
rate = rate(spikey_idx,:);
locs = locs(spikey_idx,:);

%% Do real
pre = nanmean(rate(:,change-surround:change-1),2);
post = nanmean(rate(:,change+1:change+surround),2);
achange = post-pre;
[beta,pval] = spatial_autocorrelation(thing,achange,dmin,locs,1);
true_pval = pval(1);

%% Boot
nblocks = size(rate,2);
poss_times = surround+1:nblocks-surround;
perm_beta = nan(length(poss_times),1);

% Loop over permutations
for i = 1:length(poss_times)
    
    fchange = poss_times(i);
    
    % recalculate change around this time
    fpre = nanmean(rate(:,fchange-surround:fchange-1),2);
    fpost = nanmean(rate(:,fchange+1:fchange+surround),2);

    fchange = fpost-fpre;
    [fbeta,fpval] = spatial_autocorrelation(thing,fchange,dmin,locs,1);
    perm_beta(i) = fbeta(1);
    
end

if 1
    figure
    plot(poss_times,perm_beta)
    hold on
    plot([change change],ylim)
    xlim([1 nblocks])
    title(sprintf('%s beta = %1.2f, p = %1.3f',name,beta(1),true_pval))
    print(gcf,[outfolder,name],'-dpng')
    close(gcf)
end

