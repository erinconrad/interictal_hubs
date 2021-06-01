function [mp,stp] = shaded_error_bars(times,m,st,color)

if isempty(color)
    color = [0,0.4470, 0.7410];
end

%% Flip dimensions if needed
if size(times,1) > 1
    times = times';
end

if size(m,1) > 1
    m = m';
end

if size(st,1) >1
    st = st';
end

%% Plot the line
mp = plot(times,m,'color',color,'linewidth',3);
hold on

%% Plot the patch
upper = m + st;
lower = m - st;
in_between = [upper, fliplr(lower)];
x2 = [times,fliplr(times)];
stp = fill(x2, in_between,color,'linestyle','none');
alpha(stp,0.4);

end