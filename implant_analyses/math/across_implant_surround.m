function [start,pre,post,ending] = across_implant_surround(rate,cblock,surround)

first_non_nan_pre = [];
first_non_nan_post = [];
first_non_nan_start = [];
first_non_nan_end = [];

%% Find first non-nan for beginning of file, end of file, pre cblock, postcblock
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

% Beginning of file
for i = 1:first_non_nan_pre
    if sum(~isnan(rate(:,i))) > 0 % if any electrodes are not nans
        first_non_nan_start = i; % this is the last non-nan time before the revision
        break
    end
end
    
% End of file
for i = size(rate,2)-1:-1:first_non_nan_post
    if sum(~isnan(rate(:,i))) > 0 % if any electrodes are not nans
        first_non_nan_end = i; % this is the last non-nan time before the revision
        break
    end
end

% Do some sanity checks
if isempty(first_non_nan_pre) || isempty(first_non_nan_post) || isempty(first_non_nan_start) || isempty(first_non_nan_end)
    error('why')
end

if first_non_nan_pre == first_non_nan_start || first_non_nan_post == first_non_nan_end
    error('why')
end

%% Now go surround back and forth from this
pre = max(1,first_non_nan_pre-surround):first_non_nan_pre;
post = first_non_nan_post:min(first_non_nan_post+surround,size(rate,2));
start = first_non_nan_start:first_non_nan_start + surround;
ending = first_non_nan_end-surround:first_non_nan_end;

if 0
    figure
    turn_nans_white(rate)
    hold on
    
    plot([pre(1) pre(1)],ylim,'g--','linewidth',3)
    plot([pre(end) pre(end)],ylim,'g--','linewidth',3)
    plot([start(1) start(1)],ylim,'g--','linewidth',3)
    plot([start(end) start(end)],ylim,'g--','linewidth',3)
    
    plot([post(1) post(1)],ylim,'b--','linewidth',3)
    plot([post(end) post(end)],ylim,'b--','linewidth',3)
    plot([ending(1) ending(1)],ylim,'b--','linewidth',3)
    plot([ending(end) ending(end)],ylim,'b--','linewidth',3)
    
    plot([cblock cblock],ylim,'r--','linewidth',3)
    
    pause
    close(gcf)
end
end