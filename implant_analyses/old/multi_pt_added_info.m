function multi_pt_added_info(rchange,anat,results_folder,surround,pt_names)

nch = cellfun(@(x) length(x),anat);
nch = nch';

%% Correlate number of electrodes with rate change
[r,pval] = corr(nch(~isnan(rchange)),rchange(~isnan(rchange)));

figure
plot(nch,rchange,'o','markersize',15,'linewidth',2,'color',[1 1 1]);
text(nch,rchange,pt_names,'horizontalalignment','center',...
    'fontsize',20)
xlim([min(nch) max(nch)])
ylim([min(rchange) max(rchange)])
%set(h, 'MarkerFaceColor', get(h,'Color'));
xlabel('Number of added electrodes')
ylabel('Relative change in spike rate')
title(sprintf('r = %1.2f, p = %1.3f',r,pval))
set(gca,'fontsize',20)

outfolder = [results_folder,'group_level/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end

print(gcf,[outfolder,'num_elecs'],'-dpng')

%% Correlate number of electrodes in different anatomical locations with rate change
% Group anatomical locations
loc = cell(length(anat),1);
lat = cell(length(anat),1);
all_lat = {};
all_loc = {};
for i = 1:length(anat)
    [la,lo] = anatomy_grouper(anat{i});
    loc{i} = lo;
    lat{i} = la;
    
    all_lat = [all_lat;la];
    all_loc = [all_loc;lo];
end

lat_groups = unique(all_lat);
loc_groups = unique(all_loc);

lat_n = zeros(length(lat_groups),length(anat));
loc_n = zeros(length(loc_groups),length(anat));

for i = 1:length(anat)
    for j = 1:length(lat_groups)
        % how many electrodes for patient i belong to lat group j?
        lat_n(j,i) = sum(ismember(lat{i},lat_groups{j}));
    end
    
    for j = 1:length(loc_groups)
        % how many electrodes for patient i belong to loc group j?
        loc_n(j,i) = sum(ismember(loc{i},loc_groups{j}));
    end
end

%% Find missing info
missing = loc_n(1,:) == sum(loc_n,1);
loc_n(:,missing) = [];
rchange(missing) = [];
pt_names(missing) = [];

%% Do locs
% Ignore unspecified in locs
loc_n(1,:) = [];
loc_groups(1) = [];

fig = figure;
set(fig,'position',[25 540 1400 300])
t = tiledlayout(1,4,'TileSpacing','Compact');
for k = 1:4
    nexttile
    plot(loc_n(k,:),rchange,'o','markersize',15,'linewidth',2,'color',[1 1 1]);
    text(loc_n(k,:),rchange,pt_names,'horizontalalignment','center',...
        'fontsize',20)
    [r,pval] = corr(loc_n(k,(~isnan(rchange)))',rchange(~isnan(rchange)));
    title(sprintf('%s\nr = %1.2f, p = %1.3f',loc_groups{k},r,pval))
    set(gca,'fontsize',20)
end

ylabel(t,{'Relative change','in spike rate'},'fontsize',20);
xlabel(t,{'Number of added electrodes in given location'},'fontsize',20);
print(fig,[outfolder,'anatomy',sprintf('_%d',surround)],'-dpng')


end