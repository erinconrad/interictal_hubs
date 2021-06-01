function [pre,post,cosi] = new_rate_increase(all_rate,change_block,cos,surround)

if isempty(surround)
    pre_blocks = 1:change_block-1;
    post_blocks = change_block+1:size(all_rate,2);
    cos_blocks = 1:size(cos,2);
    
else
    pre_blocks = max(1,change_block-surround):change_block-1;
    post_blocks = change_block+1:min(change_block+surround,size(all_rate,2));
    cos_blocks = 1:surround;
end

pre = nanmean(all_rate(:,pre_blocks),2);
post = nanmean(all_rate(:,post_blocks),2);

% cosi is the number of times an unchanged channel has a cospike with an
% added channel divided by the total number of spikes for that channel
cosi = nanmean(cos(:,cos_blocks)./all_rate(:,post_blocks),2);

end