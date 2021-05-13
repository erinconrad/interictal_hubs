function [spikey_idx,mean_rate] = find_spikey_elecs(all_rate,min_num,change_block,surround)

pre_blocks = max(1,change_block-surround):change_block-1;
post_blocks = change_block+1:min(change_block+surround,size(all_rate,2));

mean_rate = nanmean(all_rate(:,[pre_blocks,post_blocks]),2);
spikey_idx = mean_rate >= min_num;

end