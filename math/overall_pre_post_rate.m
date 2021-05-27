function pval = overall_pre_post_rate(rate,change)

pre_count = rate(:,1:change-1);
post_count = rate(:,change+1:end);

pre_count = nansum(pre_count,1)';
post_count = nansum(post_count,1)';

pre_count = nansum(pre_count);
post_count = nansum(post_count);



n_pre = change-1;
n_post = size(rate,2) - (change + 1) + 1;

pre_rate = pre_count/n_pre;
post_rate = post_count/n_post;

z_pt = (post_rate - pre_rate)/sqrt(pre_rate/n_pre + post_rate/n_post);
p_two = normcdf(z_pt);


end