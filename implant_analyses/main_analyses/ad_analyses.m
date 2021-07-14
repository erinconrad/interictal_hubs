function ad_analyses(whichPts,saved_out,out)


%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
main_surround = 3;
ex = 10;

%% Other info
n_surrounds = length(all_surrounds);

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

%% initialize figure
figure


set(gcf,'position',[100 100 800 650])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')



%% Power spectrum of AD
nexttile

% Find longest ad
longest = 0;
for i = 1:length(out)
    ad = out(i).ad;
    if size(ad,2) > longest
        longest = size(ad,2);
    end
end

% Get all ads
all_P = nan(length(out),ceil(longest/2));
for i = 1:length(out)
    ad = out(i).ad;
    X = nanmean(ad,1); % average across electrodes
    X(isnan(X)) = nanmedian(X); % set nans to median
    X = X-nanmean(X); % substract dc component
    X(end+1:end+longest-length(X)) = 0; % set things beyond duration of run to zero
    bdur = out(i).block_dur; % 0.5 hours
    Y = fft(X); % fft
    P = (abs(Y)).^2; % power
    fs = 1/bdur; % a freq of 2/hour
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);
    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));
    
    all_P(i,:) = P/sum(P); % normalize by total power
    if i >1
        if ~isequal(freqs,curr_freqs) % all pts should have same freqs
            error('oh no');
        end
    end
    
    curr_freqs = freqs;
    
    
end

% Turn into period and restrict to less than 100 hours
period = 1./curr_freqs;
low_period_idx = period<100;
period = period(:,low_period_idx);
all_P = all_P(:,low_period_idx);

% Do plot
m = mean(all_P,1);
[mp,stp] = shaded_error_bars(period,m,std(all_P,[],1),[]);

% Identify max
[~,I] = max(m);
plot(period(I),m(I),'*','color',[0.8500, 0.3250, 0.0980],'markersize',10,...
    'linewidth',3)
text(period(I)+2,m(I),sprintf('%1.1f hours',period(I)),'fontsize',20,...
    'horizontalalignment','left','verticalalignment','bottom');
xlabel('Period (hours)');
ylabel({'Alpha-delta ratio (AD)','normalized power spectrum'})
xlim([0 max(period)])
set(gca,'fontsize',20)
%pause(0.3)

%% Single pt example spike rate vs AD
nexttile
ad = nanmean(out(ex).ad,1);
cblock = out(ex).change_block;
sr = nansum(out(ex).rate,1)/out(ex).run_dur; % sum of all spikes across electrodes in the five minute detection period/divided by 5 minutes to get rate per minute
bdur = out(ex).block_dur;
times = (1:length(sr))*bdur;
nan_blocks = find(isnan(nanmean(out(ex).rate,1)));
ylabels = {'AD','Spikes/min'};
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980];
h = stackedplot(times,[ad;sr]','DisplayLabels',ylabels);
ax = findobj(h.NodeChildren, 'Type','Axes');
pause(0.3) % delete at your own risk
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
h.LineWidth = 2;
for k = 1:length(h.LineProperties)
    if k <= size(myColours,1)
        h.LineProperties(k).Color = myColours(k,:);
    end
end
set(gca,'fontsize',20);
xlabel('Hours')
ax = findobj(h.NodeChildren, 'Type','Axes');
set(ax, 'YTick', []);

%% All pts correlation between AD and SR
nexttile
all_corr_info = nan(length(out),2);
for i = 1:length(out)
    ad = nanmean(out(i).ad,1);
    sr = nansum(out(i).rate,1);
    
    % set nans in rate to zero
    nan_blocks = find(isnan(nanmean(out(i).rate,1)));
    sr(nan_blocks) = nan;
    ad(nan_blocks) = nan;
    
    [r,p] = corr(ad',sr','rows','pairwise');
    
    % Get n
    chLabels = out(i).unchanged_labels;
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    chLabels(ekg) = []; 
    n = length(chLabels);
    
    [z,~,~] = fisher_transform(r,n);
    all_corr_info(i,:) = [r z];
    
end


all_zs = all_corr_info(:,2);
% two sided unpaired ttest
[~,pval,~,stats] = ttest(all_zs);
rate_t = stats.tstat;
rate_df = stats.df;

% get r back by averaging the z's and z-to-r transforming. Do weighted
% average by sample size
%r = tanh(nansum(all_corr_info(:,2).*all_corr_info(:,4))./nansum(all_corr_info(:,4)));

% Plot
plot(all_corr_info(:,1),'o','linewidth',2,'markersize',10) % plot r's
hold on
ylim([-1 1])
xlim([0.5 10.5])
plot(xlim,[0 0],'k--','linewidth',2)
ylabel('AD-spike correlation')
xlabel('Patient')
xticklabels([])
xl = xlim;
yl = ylim;
if pval < 0.001
    ptext = 'p < 0.001';
else
    ptext = sprintf('p = %1.3f',pval);
end
text(xl(2),yl(2),sprintf('Mean r = % 1.2f\n%s',mean(all_corr_info(:,1)),ptext),...
    'horizontalalignment','right',...
    'verticalalignment','top','fontsize',20)
set(gca,'fontsize',20)

%
%% AD pre- and post-revision (stats for next graph)
all_p = nan(n_surrounds,1);
all_df = nan(n_surrounds,1);
all_t = nan(n_surrounds,1);
for s = 1:n_surrounds
    surround = all_surrounds(s);
    all_ad = nan(length(whichPts),2);
    for i = 1:length(whichPts)
        ad = nanmean(out(i).ad,1);
        rate = out(i).rate;
        cblock = out(i).change_block;
        [pre,post] = get_surround_times(rate,cblock,surround);

        ad_pre = nanmean(ad(pre));
        ad_post = nanmean(ad(post));
        all_ad(i,:) = [ad_pre ad_post];
    end

    % Paired ttest
    [~,pval,~,stats] = ttest(all_ad(:,1),all_ad(:,2));
    if s == main_surround
        pre_main = all_ad(:,1);
        post_main = all_ad(:,2);
    end
    all_p(s) = pval;
    all_df(s) = stats.df;
    all_t(s) = stats.tstat;
end

pval = all_p(main_surround); % which one to use for plot

%% Save table of other p-values
adT = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.2f, %s',x,y,pretty_p_text(z)),all_df,all_t,all_p,...
    'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));
writetable(adT,[main_spike_results,'ad.csv'],'WriteRowNames',true)  
%{
plot(ones(npts,1)+0.05*rand(npts,1),all_ad(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_ad(:,2),'o')
xlim([0.5 2.5])
%}

%% Also add to main supplemental table
writetable(adT,[main_spike_results,'Supplemental Table 1.xlsx'],'Range','G2:G12','WriteVariableNames',false)


%
%% AD over time all patients on graph, showing implant time in middle

nexttile
% Get full surround time
max_time_before = 0;
max_time_after = 0;
for i = 1:length(whichPts)
    ntotal = size(out(i).rate,2);
    cblock = out(i).change_block;
    nbefore = cblock-1;
    nafter = ntotal-cblock;
    if nbefore > max_time_before
        max_time_before = nbefore;
    end
    if nafter > max_time_after
        max_time_after = nafter;
    end
end

all_ad = nan(length(whichPts),max_time_before+max_time_after+1);
mid_pos = max_time_before+1;
for i = 1:length(whichPts)
    ad = nanmean(out(i).ad,1);
    cblock = out(i).change_block;
    % Put it at the correct position in the larger array. If it's the one
    % with the latest cblock, then we can start filling it up at the first
    % position. If it has an earlier cblock, then we need to pad the
    % beginning
    pos_off = mid_pos - cblock;
    all_ad(i,1+pos_off:length(ad)+pos_off) = ad;
end
bdur = out(1).block_dur;
% Combine across patients, only including times with at least 2 non nans
times = -max_time_before:max_time_after;
two_non_nans = sum(~isnan(all_ad),1)>=2;
m = nanmean(all_ad,1);
st = nanstd(all_ad,[],1);
times(~two_non_nans) = [];
m(~two_non_nans) = [];
st(~two_non_nans) = [];
shaded_error_bars(times*bdur,m,st,[]);
hold on
yl = ylim;
plot([0 0],[yl(1) 0.85],'r--','linewidth',3)
xlim([-24*8 24*8])
ylabel('AD')
xlabel('Hours relative to revision')
plot([-all_surrounds(main_surround)/2 all_surrounds(main_surround)/2],[0.87 0.87],'k-','linewidth',2)
if pval < 0.05
    text(0,0.94,sprintf('p = %1.3f',pval),'horizontalalignment','center');
else
    text(0,0.94,'ns','horizontalalignment','center','fontsize',20);
end
set(gca,'fontsize',20)


annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.9 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.42 0.1 0.1],'String','D','fontsize',30,'linestyle','none')


print(gcf,[main_spike_results,'AD'],'-dpng')
print(gcf,[main_spike_results,'AD'],'-depsc')

%% text
fprintf(['\n\nThe average correlation between the'...
    ' alpha-delta ratio and spike rate was rho = %1.2f, which was significantly '...
    'less than zero (t(%d) = %1.1f, %s).\n\n'],mean(all_corr_info(:,1)),rate_df,rate_t,ptext);

fprintf(['\n\nAcross patients, there was no difference in the alpha-delta '...
    'ratio in the pre- and post-revision periods (pre-revision mean (SD) = '...
    '%1.2f (%1.2f), post-revision mean (SD) = %1.2f (%1.2f), t(%d) = %1.1f,'...
    ' p = %1.2f).\n'],mean(pre_main),std(pre_main),mean(post_main),std(post_main),...
    all_df(main_surround),all_t(main_surround),all_p(main_surround));

end