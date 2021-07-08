function [pval,true_rho,sorted_perm_rho] = compare_vecs(rate,change,surround,nb,buffer,do_buffer)

nblocks = size(rate,2);

%% Compare spike rate pre and post change
% Get surround times, starting with first non nan
[pre,post] = get_surround_times(rate,change,surround);

pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);

vec_diff = pre_mean-post_mean;
vec_diff(isnan(vec_diff)) = [];
true_rho = vecnorm(vec_diff);

perm_rho = nan(nb,1);

for ib = 1:nb
    
    while 1
        % Make a fake change time
        fchange = randi([surround+1,nblocks-surround]);

        if do_buffer
            [pre,post] = get_surround_times(rate,fchange,surround + buffer);
        else
            [pre,post] = get_surround_times(rate,fchange,surround);
        end

        if length(pre) == 1 || length(post) == 1
            continue % bad time, try another
        end

        % recalculate
        fpre = rate(:,pre);
        fpost = rate(:,post);

        vec_diff = nanmean(fpre,2)-nanmean(fpost,2);
        vec_diff(isnan(vec_diff)) = [];
        
        perm_rho_temp = vecnorm(vec_diff);
        if isnan(perm_rho_temp)
            continue % bad time, try another
        else
            break % if not a nan, cool, accept it
        end
    end
    perm_rho(ib) = perm_rho_temp;
    
end

%% sort the permutation results
[sorted_perm_rho] = sort(perm_rho);

% this is just a one-tailed test which I think is fair
num_as_sig = sum(sorted_perm_rho>=true_rho);
pval = (num_as_sig+1)/(nb+1);

if isnan(true_rho)
    error('why');
end

if sum(isnan(perm_rho)) > 0, error('why'); end

if 0
    figure
    plot(sorted_perm_rho,'o')
    hold on
    plot(xlim,[true_rho true_rho])
    pause
    close(gcf)
end



end