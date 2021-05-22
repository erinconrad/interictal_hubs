function cosa = cosi_analysis(cos,rate,change,surround,unchanged_labels,name,...
    added_rate)

%% Get the unchanged channels that often co-spike with added channels
% The channels that, when the added channels spike, these are often also
% spiking
% Mean across blocks
cos_mean = nanmean(cos,3);

%% Find spikey added chs
added_mean = nanmean(added_rate,2);
spikey_added = added_mean > 0.5;

% Restrict cos_mean to spikey added chs
cos_mean = cos_mean(spikey_added,:);
cosa = max(cos_mean,[],1)';



%% Spikey high cos


%{
%% Get spikey chs
ch_mean = nanmean(rate,2);
spikey = ch_mean > 0.5;
cospikey_idx = find(spikey & high_cos);
high_cos = (any(cos_mean > 0.5,1))';
%% plot rate over time for these
figure
set(gcf,'position',[200 200 1000 400])
plot(rate(cospikey_idx,:)','k')
hold on
plot([change change],ylim,'r--','linewidth',2)
title(name)
pause
close(gcf)
%}

end