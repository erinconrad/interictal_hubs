function ad_analyses(whichPts,saved_out)


%% Parameters
do_test = 0;
surround = 48;
prt = [10 50];
ex = 10;
do_additional = 0;

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

%% initialize figure
figure

if do_additional
    set(gcf,'position',[100 100 1200 600])
    tiledlayout(2,3,'TileSpacing','compact','Padding','tight')
else
    set(gcf,'position',[100 100 800 650])
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
end


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
    X = nanmean(ad,1);
    X(isnan(X)) = nanmedian(X);
    X = X-nanmean(X); % substract dc
    X(end+1:end+longest-length(X)) = 0;
    bdur = out(i).block_dur;
    Y = fft(X);
    P = (abs(Y)).^2;
    fs = 1/bdur;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);
    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));
    
    all_P(i,:) = P/sum(P); % normalize by total power
    if i >1
        if ~isequal(freqs,curr_freqs)
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
%{
plot(times,ad,'linewidth',2)
hold on
plot(times,sr,'linewidth',2)
%}

%% All pts correlation between AD and SR
nexttile
all_corr_info = nan(length(out),4);
for i = 1:length(out)
    ad = nanmean(out(i).ad,1);
    sr = nansum(out(i).rate,1);
    [r,p] = corr(ad',sr','rows','pairwise');
    n = sum(~isnan(ad) & ~isnan(sr));
    [z,z_score,pval2] = fisher_transform(r,n);
    if abs(p-pval2) > 0.03
        error('oh nos');
    end
    all_corr_info(i,:) = [r z z_score n];
    
end

% Stouffer's Z-score
if sum(sum(isnan(all_corr_info)))>0, error('what'); end
z = sum(all_corr_info(:,3))/sqrt(length(out));
% get two sided p-value from the combined z score
pval = 2*normcdf(-abs(z));

% get r back by averaging the z's and z-to-r transforming. Do weighted
% average by sample size
r = tanh(nansum(all_corr_info(:,2).*all_corr_info(:,4))./nansum(all_corr_info(:,4)));

% Plot
plot(all_corr_info(:,1),'o','linewidth',2,'markersize',10)
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
text(xl(2),yl(2),sprintf('Combined r = % 1.2f\n%s',r,ptext),...
    'horizontalalignment','right',...
    'verticalalignment','top','fontsize',20)
set(gca,'fontsize',20)

%
%% AD pre- and post-revision (stats for next graph)
npts = length(whichPts);
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
%{
plot(ones(npts,1)+0.05*rand(npts,1),all_ad(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_ad(:,2),'o')
xlim([0.5 2.5])
%}

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
all_ps = nan(length(whichPts),1);
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
plot([-24 24],[0.87 0.87],'k-','linewidth',2)
if pval < 0.05
    text(0,0.94,sprintf('p = %1.3f',pval),'horizontalalignment','center');
else
    text(0,0.94,'ns','horizontalalignment','center','fontsize',20);
end
set(gca,'fontsize',20)

%% High AD - to - low AD, single patient
if do_additional 
    nexttile


    % Get difference in rate between pre and post revision
    rate = out(ex).rate;
    cblock = out(ex).change_block;
    [pre,post] = get_surround_times(rate,cblock,surround);
    rate_pre = nanmean(rate(:,pre),2);
    rate_post = nanmean(rate(:,post),2);
    re_diff_rate = rate_post - rate_pre;

    % Get high and low ad periods (remove periods in pre and post)
    ad = nanmean(out(ex).ad,1);
    idx = 1:size(ad,2);
    Y = prctile(ad,prt); % get 10th and 90th %ile AD values
    high_ad_idx = ad > Y(2) & ~ismember(idx,pre) & ~ismember(idx,post);
    low_ad_idx = ad < Y(1) & ~ismember(idx,pre) & ~ismember(idx,post);


    %{
    plot(ad)
    hold on
    plot(find(high_ad_idx),ad(high_ad_idx),'ro')
    plot(find(low_ad_idx),ad(low_ad_idx),'go')
    %}



    % Get difference in rate between high and low ad
    rate_high_ad = nanmean(rate(:,high_ad_idx),2);
    rate_low_ad = nanmean(rate(:,low_ad_idx),2);
    ad_diff_rate = rate_low_ad-rate_high_ad;
    mean_rate = nanmean(rate,2);


    % Throw out ekg and such
    unchanged_labels = out(ex).unchanged_labels;
    ekg = identify_ekg_scalp(unchanged_labels);
    ad_diff_rate = ad_diff_rate(~ekg);
    re_diff_rate = re_diff_rate(~ekg);
    mean_rate = mean_rate(~ekg);
    unchanged_labels = unchanged_labels(~ekg);

    % Corr
    if do_test
        [r,p] = corr(mean_rate,re_diff_rate,'rows','pairwise');
        plot(mean_rate,re_diff_rate,'o','color',[1 1 1])
        text(mean_rate,re_diff_rate,unchanged_labels,'horizontalalignment','center')
        xlabel('Mean rate')
        ylabel('Post-pre rate change')
    else
        [r,p] = corr(ad_diff_rate,re_diff_rate,'rows','pairwise');
        plot(ad_diff_rate,re_diff_rate,'o','color',[1 1 1])
        text(ad_diff_rate,re_diff_rate,unchanged_labels,'horizontalalignment','center')
        xlabel('Low AD-High AD rate change')
        ylabel('Post-pre rate change')
    end
    xl = xlim;
    yl = ylim;
    text(xl(2),yl(2),sprintf('r = %1.2f\np = %1.3f',r,p),...
        'horizontalalignment','right','VerticalAlignment','top')

    %% Corr all patients
    nexttile
    all_corr_info = nan(length(out),5);
    npts = length(out);
    for i = 1:npts

        % Get difference in rate between pre and post revision
        rate = out(i).rate;
        cblock = out(i).change_block;
        [pre,post] = get_surround_times(rate,cblock,surround);
        rate_pre = nanmean(rate(:,pre),2);
        rate_post = nanmean(rate(:,post),2);
        re_diff_rate = rate_post - rate_pre;

        % mean rate
        mean_rate = nanmean(rate,2);

        % Get high and low ad periods (remove periods in pre and post)
        ad = nanmean(out(i).ad,1);
        idx = 1:size(ad,2);
        Y = prctile(ad,prt); % get 10th and 90th %ile AD values
        high_ad_idx = ad > Y(2) & ~ismember(idx,pre) & ~ismember(idx,post);
        low_ad_idx = ad < Y(1) & ~ismember(idx,pre) & ~ismember(idx,post);

        % Get difference in rate between high and low ad
        rate_high_ad = nanmean(rate(:,high_ad_idx),2);
        rate_low_ad = nanmean(rate(:,low_ad_idx),2);
        ad_diff_rate = rate_low_ad-rate_high_ad;

        % Throw out ekg and such
        ekg = identify_ekg_scalp(out(i).unchanged_labels);
        ad_diff_rate = ad_diff_rate(~ekg);
        re_diff_rate = re_diff_rate(~ekg);
        mean_rate = mean_rate(~ekg);

        % Corr
        if do_test
            [r,p] = corr(re_diff_rate,mean_rate,'rows','pairwise');
            n = sum(~isnan(re_diff_rate) & ~isnan(mean_rate));
        else
            [r,p] = corr(ad_diff_rate,re_diff_rate,'rows','pairwise');
            n = sum(~isnan(ad_diff_rate) & ~isnan(re_diff_rate));
        end
        [z,z_score,pval2] = fisher_transform(r,n);
        if abs(p-pval2) > 0.03
            error('oh nos');
        end
        all_corr_info(i,:) = [r z z_score n p];

    end

    % Stouffer's Z-score
    if sum(sum(isnan(all_corr_info)))>0, error('what'); end
    z = sum(all_corr_info(:,3))/sqrt(length(out));
    % get two sided p-value from the combined z score
    pval = 2*normcdf(-abs(z));

    % get r back by averaging the z's and z-to-r transforming. Do weighted
    % average by sample size
    r = tanh(nansum(all_corr_info(:,2).*all_corr_info(:,4))./nansum(all_corr_info(:,4)));

    % Plot
    plot(all_corr_info(:,1),'o','linewidth',2)
    hold on
    ylim([-1 1])
    plot(xlim,[0 0],'k--')
    if do_test
        ylabel('Rate-revision rate diff correlation')
    else
        ylabel('AD rate diff-revision rate diff correlation')
    end
    xlabel('Patient')
    xticklabels([])
    text(9.5,0.9,sprintf('Combined r = % 1.2f\np = %1.3f',r,pval),'horizontalalignment','right')
end

annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.9 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.42 0.1 0.1],'String','D','fontsize',30,'linestyle','none')


print(gcf,[main_spike_results,'AD'],'-dpng')
print(gcf,[main_spike_results,'AD'],'-depsc')

end