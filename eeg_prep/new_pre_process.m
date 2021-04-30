function values = new_pre_process(values,which_chs)

% Just do a car (just average the non-skip chs)
values = values - repmat(mean(values(:,which_chs),2),1,size(values,2));

end