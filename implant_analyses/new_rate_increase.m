function increase = new_rate_increase(all_rate,change_block,surround)

pre_blocks = max(1,change_block-surround):change_block-1;
post_blocks = change_block+1:min(change_block+surround,size(all_rate,2));

increase = (nanmean(all_rate(:,post_blocks),2)-nanmean(all_rate(:,pre_blocks),2));

end