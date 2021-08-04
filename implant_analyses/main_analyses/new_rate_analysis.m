function new_rate_analysis(whichPts,saved_out,out)

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
%all_surrounds = 12*[2];
main_surround = 3; %24 hour peri-revision surround
main_metric = 1;
ex_p = 1;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'rate'};
n_metrics = length(all_metrics);

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
    
    %out = load([main_spike_results,'out.mat']);
    %out = out.out;
    
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

%% Get names
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

%% Do main analysis

% initialize arrays
all_ov = nan(n_metrics,n_surrounds,length(whichPts),2); % 2 is pre and post rate
all_added_elecs= nan(length(whichPts)); 
main_stats = nan(n_metrics,n_surrounds,3); % p, tstat, df
added_stats = nan(n_metrics,n_surrounds,2); % r, p

% Loop over metrics
for im = 1:n_metrics
    allcatlist = [];
    metric = all_metrics{im};
    
    % Loop over surround times
    for is = 1:n_surrounds

        surround = all_surrounds(is);
        
        %% Get individual patient data
        
        for i = 1:length(whichPts)
            
            switch metric
                case 'rate'
                    rate = out(i).rate;
                case 'ns'
                    
                    rate = out(i).metrics.ns;
                   
            end
            
            cblock = out(i).change_block;

            % Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            rate(ekg,:) = [];
            run_dur = out(i).run_dur;
            
            % Get surround times, starting with first non nan
            [pre,post] = get_surround_times(rate,cblock,surround);
            
            % Get the rate for all electrodes in the pre and post
            rate_pre = nanmean(rate(:,pre),2); % divide by 5 minutes to get rate per minute
            rate_post = nanmean(rate(:,post),2);
            
            % If doing rate, divide by run_dur to get spikes/minute
            if im == 1
                rate_pre = rate_pre/run_dur;
                rate_post = rate_post/run_dur;
            end
            
            % overall rate
            ov_pre = nansum(rate_pre); % sum across electrodes
            ov_post = nansum(rate_post); 

            % Fill up array with data
            all_ov(im,is,i,:) = [ov_pre ov_post];
            
            % Get info on added electrodes
            n_added = length(out(i).change(end).added);
            
            % Fill up info
            all_added_elecs(i,1) = n_added;
            
        end
        
        %% Do stats across patients
        
        % Does spike rate change pre-to-post?
        curr_ov = squeeze(all_ov(im,is,:,:));
        [~,p_main_rate,~,stats] = ttest(curr_ov(:,1),curr_ov(:,2));
        tstat = stats.tstat;
        df = stats.df;
        
        % (2) Does relative spike rate increase correlate with number of electrodes
        ov_change = (curr_ov(:,2)-curr_ov(:,1))./(curr_ov(:,1));
        [r_ov_all,p_ov_all] = corr(ov_change,all_added_elecs(:,1));
                
        % Fill up arrays
        main_stats(im,is,:) = [p_main_rate tstat df];
        added_stats(im,is,:) = [r_ov_all,p_ov_all];
        
    end
    
end

%% Initialize figure
im = main_metric;
is = main_surround;

figure
set(gcf,'position',[100 87 800 600])
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

% colors
cols = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.4940, 0.1840, 0.5560];

%% Example patient, full spike rate
nexttile([1 2])
rate = out(ex_p).rate;
cblock = out(ex_p).change_block;
% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(out(ex_p).unchanged_labels);
rate(ekg,:) = [];
run_dur = out(ex_p).run_dur;
block_dur = out(ex_p).block_dur;
rate = rate/run_dur;
rate = nansum(rate,1); % sum across channels
times = (1:length(rate)) * block_dur;

% Get nan blocks
nan_blocks = find(isnan(nanmean(out(ex_p).rate,1)));

% plot
plot(times,rate,'k','linewidth',1)
ylabel('Spikes/min')
 xlabel('Hour')
hold on

% Get surround times, starting with first non nan
[pre,post] = get_surround_times(out(ex_p).rate,cblock,surround);
pre = pre*block_dur;
post = post*block_dur;
xlim([0 length(rate)*block_dur]);

set(gca,'fontsize',15)
yl = ylim;
new_yl = [yl(1) 1.25*(yl(2)-yl(1))];
top = yl(1) + 1.1*(yl(2)-yl(1));
ybar = yl(1) + 1.15*(yl(2)-yl(1));
ytext=  yl(1) + 1.2*(yl(2)-yl(1));
ylim(new_yl);
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(rate,2))];
    bidx = bidx*out(ex_p).block_dur;
    
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),top,top],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
cblock = cblock*block_dur;
cp = plot([cblock cblock],[yl(1) top],'--','color',cols(3,:),'linewidth',5);
plot([pre(1) pre(end)],[ybar ybar],'color',cols(1,:),'linewidth',2);
plot([post(1) post(end)],[ybar ybar],'color',cols(2,:),'linewidth',2);

xl = xlim;
cblock_fig_units = axescoord2figurecoord(cblock,nan);
ytext_fig_units = axescoord2figurecoord(ytext,nan);
annotation('textarrow',[0.5 cblock_fig_units],...
    [ytext_fig_units 0.7],'String','Revision','color',cols(3,:),...
    'fontsize',15);

text((pre(1)+pre(end))/2,ytext,'Pre','color',cols(1,:),'fontsize',15,...
    'HorizontalAlignment','center')
text((post(1)+post(end))/2,ytext,'Post','color',cols(2,:),'fontsize',15,...
    'HorizontalAlignment','center')
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('Patient %d',ex_p),'fontsize',15,'VerticalAlignment','Top')


%% Spike rate change across electrodes
nexttile
curr_ov = squeeze(all_ov(im,is,:,:));
plot([1 2],curr_ov','linewidth',2,'color',cols(1,:))
hold on
xlim([0 3])
xticks([1 2])
xticklabels({'Pre','Post'})
ylabel('Spikes/min')

% Add stats
yl = ylim;
new_yl = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
ybar = yl(1) + 1.1*(yl(2)-yl(1));
yp = yl(1) + 1.15*(yl(2)-yl(1));
ylim(new_yl);
plot([1 2],[ybar ybar],'k')
text(1.5,yp,sprintf('%s',get_asterisks(main_stats(im,is,1),1)),...
    'horizontalalignment','center','fontsize',15);
set(gca,'fontsize',15)

%% Spike rate correlated with number of added electrodes
nexttile
curr_ov = squeeze(all_ov(im,is,:,:));
ov_change = (curr_ov(:,2)-curr_ov(:,1))./(curr_ov(:,1));
plot(all_added_elecs(:,1),ov_change,'o','markersize',15,'linewidth',2)
xlabel('# Added electrodes');
ylabel('Relative spike rate change')
set(gca,'fontsize',15)
ylim([-1 2.5])
xl = xlim;
yl = ylim;
text(xl(2),yl(1),sprintf('r = %1.2f, %s',...
    added_stats(im,is,1),get_p_text(added_stats(im,is,2))),...
    'horizontalalignment','right','fontsize',15)



end