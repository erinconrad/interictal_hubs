function rate_all_pts(whichPts,saved_out)

%% Parameters
surround = 48;
nb = 1e4;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

if saved_out == 1
    
    out = load([main_spike_results,'out.mat']);
    out = out.out;
    
else
    out = initialize_out_struct(length(whichPts));
    
    %% Get spike details
    fprintf('Getting spike details for pt...\n');
    for i = 1:length(whichPts)
        p = whichPts(i);
        fprintf('%d of %d\n',i,length(whichPts));
        out(i) = get_gdf_details(p);
    end
    save([main_spike_results,'out'],'out');
end

%% Get overall rate and do significance testing
for i = 1:length(whichPts)
    out(i).overall_rate = nansum(out(i).rate,1);
    out(i).nan_blocks = find(isnan(nanmean(out(i).rate,1)));
    
    out(i).overall_rate_pval = mc_overall_rate(out(i).overall_rate,...
        surround,out(i).change_block,nb);
end

%% Initialize figure
figure
set(gcf,'position',[252 547 1001 250])
tiledlayout(1,2)

%% Example rate
nexttile
p = 1;
curr_rate = out(p).overall_rate /out(p).run_dur;
curr_times = (1:length(curr_rate)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
plot(curr_times,curr_rate,'k','linewidth',2)
hold on
nan_blocks = out(i).nan_blocks;

xlim([0 length(curr_rate)*out(p).block_dur]);
ylabel('Spikes/min')
xlabel('Hour')
title(sprintf('%Pt %d',p));
set(gca,'fontsize',20)
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(curr_rate,2))];
    bidx = bidx*out(p).block_dur;
    yl = ylim;
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),yl(2),yl(2)],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cp = plot([curr_change curr_change],ylim,'r--','linewidth',4);
legend([cp ap],{'Revision','Data missing'},'fontsize',20,'location','northeast')

%% Change for each patient
nexttile
for i = 1:length(whichPts)
    cb = out(i).change_block;
    rate = out(i).overall_rate/out(i).run_dur;
end




end