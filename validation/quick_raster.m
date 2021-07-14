function quick_raster(out,p)
%{
This is a handy validation function. You feed in the out file and specify
the patient, and it plots a raster of spike rate in different channels over
time. You click on the block and electrode of interest and it picks a
random spike from that block-electrode pair, downloads the data from ieeg,
and shows you the spike.
%}

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];
name = out(p).name;

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
chLabels = out(p).unchanged_labels;
cblock = out(p).change_block;

ekg = identify_ekg_scalp(chLabels);

rate = rate(~ekg,:);
chLabels = chLabels(~ekg);
findices = out(p).findices;
bindices = out(p).bindices;


figure
turn_nans_white(rate)
hold on
plot([cblock cblock],ylim,'r--','linewidth',4)
yticks(1:length(chLabels))
yticklabels(chLabels)

%
while 1
    try
        [x,y] = ginput;
    catch
        break
    end
    chLab = chLabels{round(y(end))};
    fidx = findices(round(x(end)));
    bidx = bindices(round(x(end)));
    fprintf('\nShowing spikes for %s ch %s file %d block %d\n',...
        pt(pt_p).name,chLab,fidx,bidx);

    plot_spikes_by_ch(pt_p,chLab,fidx,bidx,spikes)

end
%}


end