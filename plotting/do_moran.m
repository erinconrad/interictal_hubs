function do_moran(rate_inc,unchanged_locs,spikey_idx,added_locs,name,...
    results_folder,run_dur,rate,change,surround)

only_spikey = 1;

%% Sig
pval = moran_perm_analysis(rate,change,surround,unchanged_locs,spikey_idx);

figure
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),100,'k');
hold on
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),100,'p','r');

if only_spikey
    spikey_locs = unchanged_locs(spikey_idx,:);
    spikey_rate_inc = rate_inc(spikey_idx);
    rate = rate(spikey_idx,:);
else
    spikey_locs = unchanged_locs;
    spikey_rate_inc = rate_inc;
end


scatter3(spikey_locs(:,1),spikey_locs(:,2),spikey_locs(:,3),100,spikey_rate_inc/run_dur,...
'filled');

%% Moran
dmin = ceil(nanmedian(vecnorm(diff(unchanged_locs,1),2,2)));
MI = moran_index(spikey_locs,spikey_rate_inc,dmin);
title(sprintf('%s MI = %1.2f, p = %1.3f, comp p = %1.3f',name,MI.I,MI.p,pval))
h = colorbar;
ylabel(h,'Rate change (spikes/min)');
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'fontsize',20)

%% Output folder
outfolder = [results_folder,'moran/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
print(gcf,[outfolder,name,'_moran'],'-dpng')
close(gcf)

end