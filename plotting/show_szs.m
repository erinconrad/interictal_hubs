function show_szs(out,pt,p)

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'job_talk/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder);
end

%% Get spike rate out
rate = out(p).rate./out(p).run_dur;
curr_times = (1:size(rate,2)) * out(p).block_dur/24;
ekg = identify_ekg_scalp(out(p).unchanged_labels);
rate(ekg,:) = [];

%% Remove all nan elecs
rate(sum(~isnan(rate),2)==0,:) = [];

%% Get corresponding pt index
found = 0;
for j = 1:length(pt)
    if strcmp(out(p).name,pt(j).name)
        found = 1;
        break
    end
end
if found == 0, error('oh no'); end

%% Get sz info
all_szs = [];
nblocks_prior = 0;

for f = 1:length(pt(j).ieeg.file)
    block_times = mean(pt(j).ieeg.file(f).block_times,2);
    for s = 1:length(pt(j).ieeg.file(f).sz)
        sz_time = pt(j).ieeg.file(f).sz(s).ueo;
        % find closest block
        [~,closest_block] = min(abs(sz_time-block_times));
        
        all_szs = [all_szs;closest_block + nblocks_prior];
    end
    nblocks_prior = nblocks_prior + length(block_times);
    if f == length(pt(j).ieeg.file)
        if nblocks_prior ~= size(rate,2), error('what');end
    end
end
all_szs = all_szs* out(p).block_dur/24;

%% plot
figure
set(gcf,'position',[100 87 1000 500])
h = turn_nans_white(rate);
set(h,'XData',[0:curr_times(end)]);
xlim([0 10])
hold on
for i = 1:length(all_szs)
    plot([all_szs(i) all_szs(i)],ylim,'r--','linewidth',2);
end
xlabel('Day')
ylabel('Electrode')
yticklabels([])
set(gca,'fontsize',20)
c = colorbar;
ylabel(c,'Spikes/elec/min','fontsize',20);

print([out_folder,'raster'],'-dpng')

end