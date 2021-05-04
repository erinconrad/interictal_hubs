function [rate_increase,spikey_labels,idx_above_min] = find_rate_increase(all_rate,change_block,labels)

%% Parameters
min_rate = 1;
pre_blocks = 1:change_block-1;
post_blocks = change_block+1:min(change_block+100,size(all_rate,2));

%% Get indices of electrodes above the mean rate
mean_rate = nanmean(all_rate,2);
idx_above_min = mean_rate > min_rate;
spikey_labels = labels(idx_above_min);
spikey_rate = all_rate(idx_above_min,:);

%% Get relative spike increase pre-to-post
rate_increase = (nanmean(spikey_rate(:,post_blocks),2)-nanmean(spikey_rate(:,pre_blocks),2))./nanmean(spikey_rate(:,pre_blocks),2);

%% Sort by biggest rate increase
[big_rate_increase,I] = sort(rate_increase,'descend');

T = table(spikey_labels(I),big_rate_increase);


end