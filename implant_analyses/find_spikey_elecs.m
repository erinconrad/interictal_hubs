function [good_idx,mean_rate] = find_spikey_elecs(all_rate,min_num,change_block,surround)

pre_blocks = max(1,change_block-surround):change_block-1;
post_blocks = change_block+1:min(change_block+surround,size(all_rate,2));

mean_rate = nanmean(all_rate(:,[pre_blocks,post_blocks]),2);
spikey_idx = mean_rate >= min_num;

% Also remove if lots of nans
perc_nans = sum(isnan(all_rate(:,[pre_blocks,post_blocks])),2)/length([pre_blocks,post_blocks]);
nany = perc_nans > 0.5;

good_idx = spikey_idx & ~nany;

end