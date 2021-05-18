function pval = ec_permanova(rate,spikey_idx,change,surround,relative)

rate(~spikey_idx,:) = [];
nblocks = size(rate,2);
nb = 1e3;


%% Is the pre-post difference bigger than other changes over time?
if relative

    %% Compute difference pre- and post
    pre = rate(:,change-surround:change-1);
    post = rate(:,change+1:change+surround);

    % Mean across pre and post time periods
    a = nanmean(post,2) - nanmean(pre,2);

    % vecnorm of this difference
    true_diff = vecnorm(a(~isnan(a)));
    perm_diff = nan(nb,1);
    for ib = 1:nb

        % Make a fake change time
        fchange = randi([surround+1,nblocks-surround]);

        % recalculate
        fpre = rate(:,fchange-surround:fchange-1);
        fpost = rate(:,fchange+1:fchange+surround);

        a = nanmean(fpost,2) - nanmean(fpre,2);
        perm_diff(ib) = vecnorm(a(~isnan(a)));

    end

    %% sort the permutation results
    [sorted_perm_diff] = sort(perm_diff);

    if 0
        figure
        plot(sorted_perm_diff,'o')
        hold on
        plot(xlim,[true_diff true_diff])

    end

    %% sort the permutation results
    [sorted_perm_diff] = sort(perm_diff);
    num_as_sig = sum(sorted_perm_diff>=true_diff);
    pval = (num_as_sig+1)/(nb+1);
    
else
    %% Is there a pre-post difference at all? (Simple permanova)
    %% Compute difference pre- and post
    pre = rate(:,change-surround:change-1);
    post = rate(:,change+1:change+surround);
    
    %% concatenate, keep indices saying what's pre and what's post
    all = [pre,post];
    idx = [ones(size(pre,2),1);zeros(size(post,2),1)];
    
    %% Get dissimilarity matrix
    dis = get_dis(all);
    
    result = f_permanova(dis,idx,nb,1);
    pval = result.p;
    
end

end

function dis = get_dis(X)

dis = nan(size(X,2),size(X,2));

% Loop over columns
for i = 1:size(X,2)
    for j = 1:size(X,2)
        a = X(:,i)-X(:,j);
        dis(i,j) = vecnorm(a(~isnan(a)));
    end
end


end