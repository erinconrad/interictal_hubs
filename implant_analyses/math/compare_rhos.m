function [pval,true_rho,perm_rho,n] = compare_rhos(rate,change,surround,nb,type,buffer,do_buffer,real_rate)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
% Get surround times, starting with first non nan
[pre,post] = get_surround_times(real_rate,change,surround);
pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);

% Get spike rate stability (Spearman correlation pre-to-post)
true_rho = corr(pre_mean,post_mean,'Type',type,'rows','complete');
n = sum(~isnan(pre_mean) & ~isnan(post_mean));

perm_rho = nan(nb,1);

for ib = 1:nb
    
    while 1
    
        % Make a random pseudo-revision time
        fchange = randi([surround+1,nblocks-surround]);
       
        if do_buffer
            [pre,post] = get_surround_times(real_rate,fchange,surround+buffer);
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
        perm_rho_temp = corr(nanmean(fpre,2),nanmean(fpost,2),'Type',type,'rows','complete');
        % pairwise and complete seem to give same result (i believe those
        % are only different if the input is a matrix with more than 2
        % columns)
        
        if isnan(perm_rho_temp)
            continue % bad time, try another
        else
            break
        end
        
    end
    
    perm_rho(ib) = perm_rho_temp;
    
end

%% sort the permutation results
[sorted_perm_rho] = sort(perm_rho);

% this is just a one-tailed test (because a two-tailed test doesn't make
% sense)
num_as_sig = sum(sorted_perm_rho<=true_rho); % number of mc corr's less than or equal to true corr
pval = (num_as_sig+1)/(nb+1);


% I should not have any nans
if sum(isnan(perm_rho))>0
    error('why');
end

if isnan(true_rho)
    error('why')
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