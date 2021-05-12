function [rate_increase,spikey_labels,idx_above_min,mean_rate_post,...
    big_inc_labels,big_inc_rate_sorted,abs_increase] = find_rate_increase(all_rate,change_block,labels)

%% Parameters
do_median = 1;
min_rate = 0.5;
big_inc_thresh = 1;
num_blocks = 100;
pre_blocks = max(1,change_block-num_blocks):change_block-1;%1:change_block-1;
post_blocks = change_block+1:min(change_block+num_blocks,size(all_rate,2));
%pre_blocks = 1:change_block-1;
%post_blocks = change_block:size(all_rate,2);

%% Get indices of electrodes above the mean rate
if do_median
    mean_rate = nanmedian(all_rate(:,post_blocks),2);
else
    mean_rate = nanmean(all_rate(:,post_blocks),2);
end

idx_above_min = mean_rate > min_rate;
spikey_labels = labels(idx_above_min);
spikey_rate = all_rate(idx_above_min,:);

%% Get relative spike increase pre-to-post
if do_median
    rate_increase = (nanmedian(spikey_rate(:,post_blocks),2)-nanmedian(spikey_rate(:,pre_blocks),2))./nanmedian(spikey_rate(:,pre_blocks),2);
    mean_rate_post = nanmedian(spikey_rate(:,post_blocks),2);
    abs_increase = (nanmedian(spikey_rate(:,post_blocks),2)-nanmedian(spikey_rate(:,pre_blocks),2));
else
    rate_increase = (nanmean(spikey_rate(:,post_blocks),2)-nanmean(spikey_rate(:,pre_blocks),2))./nanmean(spikey_rate(:,pre_blocks),2);
    mean_rate_post = nanmean(spikey_rate(:,post_blocks),2);
    abs_increase = (nanmean(spikey_rate(:,post_blocks),2)-nanmean(spikey_rate(:,pre_blocks),2));
end

idx = rate_increase > big_inc_thresh;
big_inc_rate = rate_increase(idx);
[big_inc_rate_sorted,I] = sort(big_inc_rate,'descend');
big_inc_labels = spikey_labels(idx);
big_inc_labels = big_inc_labels(I);



end