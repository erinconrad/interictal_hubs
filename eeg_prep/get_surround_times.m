function [pre,post] = get_surround_times(rate,cblock,surround)

first_non_nan_pre = [];
first_non_nan_post = [];


%% Find the first non nan pre and post-cblock time
% Pre
for i = cblock-1:-1:1
    if sum(~isnan(rate(:,i))) > 0 % if any electrodes are not nans
        first_non_nan_pre = i; % this is the last non-nan time before the revision
        break
    end
end

% Post
for i = cblock+1:size(rate,2)
    if sum(~isnan(rate(:,i))) > 0
        first_non_nan_post = i; % first non-nan time after revision
        break
    end
end

if isempty(first_non_nan_pre) || isempty(first_non_nan_post)
    pre = nan;
    post = nan;
    return;
end

%% Now go surround back and forth from this
pre = max(1,first_non_nan_pre-surround):first_non_nan_pre;
post = first_non_nan_post:min(first_non_nan_post+surround,size(rate,2));

if 0
    figure
    turn_nans_white(rate)
    hold on
    plot([cblock cblock],ylim,'r--','linewidth',3)
    plot([pre(1) pre(1)],ylim,'g--','linewidth',3)
    plot([pre(end) pre(end)],ylim,'g--','linewidth',3)
    
    plot([post(1) post(1)],ylim,'b--','linewidth',3)
    plot([post(end) post(end)],ylim,'b--','linewidth',3)
    
    pause
    close(gcf)
end

end