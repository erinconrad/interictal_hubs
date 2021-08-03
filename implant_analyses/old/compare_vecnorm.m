function [pval,true_dist,perm_dist,n] = compare_vecnorm(rate,change,surround,nb,buffer,do_buffer,real_rate)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
% Get surround times, starting with first non nan (no buffer)
[pre,post] = get_surround_times(real_rate,change,surround);
pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);

% Get spike rate difference (Spearman correlation pre-to-post)
vec_diff = pre_mean-post_mean;
vec_diff(isnan(vec_diff)) = 0;
true_dist = vecnorm(vec_diff);
n = sum(~isnan(pre_mean) & ~isnan(post_mean));

perm_dist = nan(nb,1);

for ib = 1:nb
    
    while 1
    
        % Make a random pseudo-revision time
        fchange = randi([surround+1,nblocks-surround]);
       
        % For the Monte Carlo simulation, make the surround somewhat longer
        % to account for the gap in recording
        if do_buffer
            [pre,post] = get_surround_times_buffer(real_rate,fchange,surround,buffer);
        else
            [pre,post] = get_surround_times(real_rate,fchange,surround);
        end

        if length(pre) == 1 || length(post) == 1
            
            continue % bad time, try another
        end

        % recalculate
        fpre = rate(:,pre);
        fpost = rate(:,post);
        
        % Calculate fake spike stability
        vec_diff = nanmean(fpre,2)-nanmean(fpost,2);
        vec_diff(isnan(vec_diff)) = 0;
        perm_dist_temp = vecnorm(vec_diff);
        % pairwise and complete seem to give same result (i believe those
        % are only different if the input is a matrix with more than 2
        % columns)
        
        if isnan(perm_dist_temp)
            continue % bad time, try another
        else
            break
        end
        
    end
    
    perm_dist(ib) = perm_dist_temp;
    
end

%% sort the permutation results
[sorted_perm_dist] = sort(perm_dist);

% this is just a one-tailed test (because a two-tailed test doesn't make
% sense)
num_as_sig = sum(sorted_perm_dist>=true_dist); % number of mc corr's less than or equal to true corr
pval = (num_as_sig+1)/(nb+1);


% I should not have any nans
if sum(isnan(perm_dist))>0
    error('why');
end

if isnan(true_dist)
    error('why')
end

if 0
    figure
    plot(sorted_perm_dist,'o')
    hold on
    plot(xlim,[true_dist true_dist])
    title(sprintf('p = %1.3f',pval))
    pause
    close(gcf)
end



end