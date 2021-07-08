function [true_cat,pseudo_cats] = mc_cat(rate,change,surround,nb)

nblocks = size(rate,2);
 % Get surround times, starting with first non nan
[pre,post] = get_surround_times(rate,change,surround);
pre_mean = nanmean(rate(:,pre),2);
post_mean = nanmean(rate(:,post),2);
nlist = length(pre_mean);

% Get true cat
true_cat = cat_plot(pre_mean,post_mean);

pseudo_cats = nan(nb,nlist);

for ib = 1:nb
    
    if mod(ib,100) == 0
        ib
    end
    
    while 1
    
        % Make a fake change time
        fchange = randi([surround+1,nblocks-surround]);
        [pre,post] = get_surround_times(rate,fchange,surround);

        if length(pre) == 1 || length(post) == 1 
            continue % bad time, try another
        end
        
        % recalculate
        fpre = rate(:,pre);
        fpost = rate(:,post);
        
        pseudo_cats_temp = cat_plot(fpre,fpost);
        
    end
    
    pseudo_cats(ib,:) = pseudo_cats_temp;
    
    
end

end