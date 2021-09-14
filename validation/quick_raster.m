function quick_raster(out,p,just_save)
%{
This is a handy validation function. You feed in the out file and specify
the patient, and it plots a raster of spike rate in different channels over
time. You click on the block and electrode of interest and it picks a
random spike from that block-electrode pair, downloads the data from ieeg,
and shows you the spike.
%}

%% parameters
do_alt = 0;
surround = 24;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
if do_alt
    spike_folder = [results_folder,'alt/spikes/'];
    out_folder = [results_folder,'alt/raster/'];
else
    spike_folder = [results_folder,'new_spikes/'];
    out_folder = [results_folder,'raster/'];
end
name = out(p).name;

if exist(out_folder,'dir')==0,mkdir(out_folder);end

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% get pt-p
pt_p = nan;
for i = 1:length(pt)
    if strcmp(pt(i).name,name)
        pt_p = i;
        break
    end
end


spikes = load([spike_folder,sprintf('%s_spikes.mat',name)]);
spikes = spikes.spikes;

if ischar(p)
    for i = 1:length(out)
        if strcmp(p,out(i).name)
            p = i;
            break
        end
    end
end
        

rate = out(p).rate;
rate_added = out(p).rate_added;
added_labels = out(p).added_labels;
chLabels = out(p).unchanged_labels;
cblock = out(p).change_block;

ekg = identify_ekg_scalp(chLabels);

rate = rate(~ekg,:);
chLabels = chLabels(~ekg);
findices = out(p).findices;
bindices = out(p).bindices;

[pre,post] = get_surround_times(rate,cblock,surround);
% Get relative rate change
pre_rate = nanmean(rate(:,pre),2);
post_rate = nanmean(rate(:,post),2);
rel_rate_change = (post_rate-pre_rate)./abs(pre_rate);

[~,I] = sort(rel_rate_change,'descend');
table(chLabels(I),rel_rate_change(I))

all_rate = [rate;rate_added];
all_labels = [chLabels;added_labels];
nunchanged = length(chLabels);

figure
set(gcf,'position',[10 10 1100 800])
turn_nans_white(all_rate)
hold on
plot([cblock cblock],ylim,'r--','linewidth',4)
plot(xlim,[nunchanged nunchanged],'r--','linewidth',4);
yticks(1:length(all_labels))
yticklabels(all_labels)
title(out(p).name)

if just_save
    print(gcf,[out_folder,name],'-dpng');
    close(gcf)
    
else
    while 1
        try
            [x,y] = ginput;
        catch
            break
        end
        chLab = all_labels{round(y(end))};
        fidx = findices(round(x(end)));
        bidx = bindices(round(x(end)));
        fprintf('\nShowing spikes for %s ch %s file %d block %d\n',...
            pt(pt_p).name,chLab,fidx,bidx);

        plot_spikes_by_ch(pt_p,chLab,fidx,bidx,spikes)

    end
end


end