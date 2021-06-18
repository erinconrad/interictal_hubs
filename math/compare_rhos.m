function pval = compare_rhos(rate,change,surround,nb,only_pre,type)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
% Get surround times, starting with first non nan
[pre,post] = get_surround_times(rate,change,surround);
%pre = rate(:,change-surround:change-1);
%post = rate(:,change+1:change+surround);
pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);
true_rho = corr(pre_mean,post_mean,'Type',type,'rows','pairwise');

perm_rho = nan(nb,1);

for ib = 1:nb
    
    % Make a fake change time
    if only_pre
        fchange = randi([surround+1,change-1]);
    else
        fchange = randi([surround+1,nblocks-surround]);
    end
    
    [pre,post] = get_surround_times(rate,fchange,surround);
    
    if length(pre) == 1 || length(post) == 1
        perm_rho(ib) = nan;
        continue
    end
    
    % recalculate
    fpre = rate(:,pre);
    fpost = rate(:,post);
    
    perm_rho(ib) = corr(nanmean(fpre,2),nanmean(fpost,2),'Type',type,'rows','pairwise');
    
end

%% sort the permutation results
[sorted_perm_rho] = sort(perm_rho);

% this is just a one-tailed test which I think is fair
num_as_sig = sum(sorted_perm_rho<=true_rho);
pval = (num_as_sig+1)/(nb+1);

if isnan(true_rho)
    pval = nan;
end

if 0
    figure
    plot(sorted_perm_rho,'o')
    hold on
    plot(xlim,[true_rho true_rho])
    pause
    close(gcf)
end



end