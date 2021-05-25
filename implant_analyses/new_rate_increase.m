function [pre,post,cosi] = new_rate_increase(all_rate,change_block,cos,surround)

pre_blocks = max(1,change_block-surround):change_block-1;
post_blocks = change_block+1:min(change_block+surround,size(all_rate,2));

pre = nanmean(all_rate(:,pre_blocks),2);
post = nanmean(all_rate(:,post_blocks),2);
cos_post = nanmean(cos(:,1:surround),2);

% cosi is the number of times an unchanged channel has a cospike with an
% added channel divided by the total number of spikes for that channel
cosi = nanmean(cos(:,1:surround)./all_rate(:,post_blocks),2);

end