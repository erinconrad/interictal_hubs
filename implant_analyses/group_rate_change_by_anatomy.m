function group_rate_change_by_anatomy(ana_lat,ana_loc,rate_increase,name,results_folder,surround)



%% Group by localization
% group groups
[ana_loc_groups,~,ic] = unique(ana_loc);
nLocs = length(ana_loc_groups);
grouped_rates_loc = cell(nLocs,1);
rates_table = [];
loc_table = {};

% loop through anatomy groups
for i = 1:nLocs
    
    % get the channels corresponding to that anatomy group
    chs = find(ic == i);
    
    % Get the average rate change for those chs
    grouped_rates_loc{i} = (rate_increase(chs));
    rates_table = [rates_table;rate_increase(chs)];
    loc_table = [loc_table;ana_loc(chs)];
end

% Prep table for anova
[p_loc,~,~] = anova1(rates_table,loc_table,'off');

%% Group by laterality
% group groups
[ana_lat_groups,~,ic] = unique(ana_lat);
nLats = length(ana_lat_groups);
grouped_rates_lat = cell(nLats,1);
rates_table = [];
loc_table = {};

% loop through anatomy groups
for i = 1:nLats
    
    % get the channels corresponding to that anatomy group
    chs = find(ic == i);
    
    % Get the average rate change for those chs
    grouped_rates_lat{i} = (rate_increase(chs));
    
    rates_table = [rates_table;rate_increase(chs)];
    loc_table = [loc_table;ana_lat(chs)];
end

% Prep table for anova
[p_lat,~,~] = anova1(rates_table,loc_table,'off');

%% Plot
figure
set(gcf,'position',[47 418 1394 380])
subplot(1,2,1)
for i = 1:length(grouped_rates_loc)
    plot(i+0.05*rand(length(grouped_rates_loc{i}),1),grouped_rates_loc{i},'o',...
        'markersize',15)
    hold on
end
xticks(1:nLocs)
xticklabels(ana_loc_groups)
xtickangle(45)
xlim([0.5 length(grouped_rates_loc)+0.5])
set(gca,'fontsize',15)
ylabel({'Median spike rate increase','after revision (spikes/min)'})
title(sprintf('%s localization ANOVA p = %1.3f',name,p_loc))

subplot(1,2,2)
for i = 1:length(grouped_rates_lat)
    plot(i+0.05*rand(length(grouped_rates_lat{i}),1),grouped_rates_lat{i},'o',...
        'markersize',15)
    hold on
end
xlim([0.5 length(grouped_rates_lat)+0.5])
xticks(1:nLats)
xticklabels(ana_lat_groups)
xtickangle(60)
set(gca,'fontsize',15)
title(sprintf('%s laterality ANOVA p = %1.3f',name,p_lat))

%% Prep output folder
out_folder = [results_folder,'anatomy/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder);
end

print(gcf,[out_folder,name,'_ana'],'-dpng');
close(gcf)


end