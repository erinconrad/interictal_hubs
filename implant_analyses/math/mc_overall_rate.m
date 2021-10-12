function pval = mc_overall_rate(rate,surround,change,nb,full_rate,buffer,do_buffer)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
[pre,post,pre_nans,post_nans] = get_surround_times(full_rate,change,surround);
pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);

pre_mean = nanmean(pre_mean);
post_mean = nanmean(post_mean);
true_diff = post_mean - pre_mean;

perm_diff = nan(nb,1);

for ib = 1:nb
    
    while 1
        % Make a fake change time
        fchange = randi([surround+1,nblocks-surround]);

        % For the Monte Carlo simulation, make the surround somewhat longer
        % to account for the gap in recording
        if do_buffer
            [fpre,fpost] = get_surround_times_buffer(full_rate,fchange,surround,buffer,pre_nans,post_nans);
        else
            [fpre,fpost] = get_surround_times(full_rate,fchange,surround);
        end

        if length(fpre) == 1 || length(fpost) == 1, continue; end
        
        fpre = nanmean(rate(:,fpre),2);
        fpost = nanmean(rate(:,fpost),2);

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

if 0
    figure
    plot(sorted_perm_diff,'o')
    hold on
    plot(xlim,[true_diff true_diff])
    title(sprintf('p = %1.3f',pval))
    pause
    close(gcf)
end


end