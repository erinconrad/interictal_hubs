function do_moran(rate_inc,unchanged_locs,spikey_idx,added_locs,name,results_folder)

only_spikey = 0;

figure
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),100,'k');
hold on
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),100,'p','r');

if only_spikey
    spikey_locs = unchanged_locs(spikey_idx,:);
    spikey_rate_inc = rate_inc(spikey_idx);
else
    spikey_locs = unchanged_locs;
    spikey_rate_inc = rate_inc;
end


scatter3(spikey_locs(:,1),spikey_locs(:,2),spikey_locs(:,3),100,spikey_rate_inc,...
'filled');

%% Moran
dmin = ceil(nanmedian(vecnorm(diff(unchanged_locs,1),2,2)));
MI = moran_index(spikey_locs,spikey_rate_inc,dmin);
title(sprintf('%s MI = %1.2f p = %1.3f',name,MI.I,MI.p))

%% Output folder
outfolder = [results_folder,'moran/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
print(gcf,[outfolder,name,'_moran'],'-dpng')
close(gcf)

end