function new_rate_analysis(whichPts,saved_out,out)

%{
This analysis looks for 1) a change in overall spike rate pre-to-post
revision, 2) whether the change in spike rate distribution across
electrodes is larger peri-revision than expected by chance, and 3) whether
the change in connectivity distribution across electrodes is larger
peri-revision than expected by chance.
%}

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
main_surround = 3; %24 hour peri-revision surround
main_metric = 1;
ex_p = 4;
nb = 1e4;%1e4; % CHANGE

%% Other info
do_buffer = 1;
n_surrounds = length(all_surrounds);
all_metrics = {'rate','ns'};
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
npts = length(whichPts);

%% Do main analysis

% initialize arrays
all_ov = nan(n_metrics,n_surrounds,length(whichPts),2); % 2 is pre and post rate
all_added_elecs= nan(length(whichPts)); 
main_stats = nan(n_metrics,n_surrounds,3); % p, tstat, df
added_stats = nan(n_metrics,n_surrounds,2); % r, p

% ROS arrays
all_all_true_r = nan(n_surrounds,length(whichPts),n_metrics);
all_all_mc_r = nan(n_surrounds,length(whichPts),nb,n_metrics);
all_all_p = nan(n_surrounds,n_metrics);
all_all_all_p = nan(n_surrounds,length(whichPts),n_metrics);
all_rate_pval = nan(n_surrounds,length(whichPts),n_metrics);

% Loop over metrics
for im = 1:n_metrics
    allcatlist = [];
    metric = all_metrics{im};
    
    % Loop over surround times
    for is = 1:n_surrounds

        surround = all_surrounds(is);
        
        % Initialize correlation stuff
        all_ps = nan(length(whichPts),1);
        
        %% Get individual patient data
        
        for i = 1:length(whichPts)
            
            switch metric
                case 'rate'
                    rate = out(i).rate;
                case 'ns'
                    
                    rate = out(i).metrics.ns;
                   
            end
            
            cblock = out(i).change_block;
            
            % Get file gap
            buffer = file_gaps(out(i).name);

            % Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            rate(ekg,:) = [];
            run_dur = out(i).run_dur;
            
            % Get surround times, starting with first non nan
            [pre,post] = get_surround_times(rate,cblock,surround);
            
            % Get the rate for all electrodes in the pre and post. Note
            % that by first taking the average across times and then across
            % electrodes, I get a slightly different result than if I take
            % the mean across times and electrodes simultaneously due to
            % the presence of nans. There aren't that many nans so it
            % doesn't make a huge difference. But by doing it this way I am
            % arguably overweighting the non-nan times on the electrodes
            % with frequent nans.
            rate_pre = nanmean(rate(:,pre),2); % divide by 5 minutes to get rate per minute
            rate_post = nanmean(rate(:,post),2);
            
            % If doing rate, divide by run_dur to get spikes/minute
            if im == 1
                rate_pre = rate_pre/run_dur;
                rate_post = rate_post/run_dur;
            end
            
            % overall rate
            ov_pre = nanmean(rate_pre); % mean across electrodes
            ov_post = nanmean(rate_post); 
            
            % Get stats on individual patient overall rate change
            rate_pval = mc_overall_rate(rate,surround,cblock,nb,out(i).rate,buffer,do_buffer);
            all_rate_pval(is,i,im) = rate_pval;
            
            % Fill up array with data
            all_ov(im,is,i,:) = [ov_pre ov_post];
            
            % Get info on added electrodes
            n_added = length(out(i).change(end).added);
            
            % Fill up info
            all_added_elecs(i,1) = n_added;
            
            % Get rho and statistics (in the math folder)
            [pval_curr,true_rho,mc_rho] = ...
                compare_rhos(rate,cblock,surround,nb,'Spearman',buffer,do_buffer,out(i).rate);

            all_all_mc_r(is,i,:,im) = mc_rho;
            all_all_true_r(is,i,im) = true_rho;
            all_all_all_p(is,i,im) = pval_curr;
            all_ps(i) = pval_curr;
            
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
        
        % Fisher test to combine pvalues for ros
        if sum(isnan(all_ps)) > 0, error('why'); end
        X_2 = -2 * sum(log(all_ps));
        sum_p = 1-chi2cdf(X_2,2*length(all_ps));

        all_all_p(is,im) = sum_p;
        
    end
    
end

%% Initialize figure
im = main_metric;
is = main_surround;

figure
set(gcf,'position',[100 1 700 600])
tiledlayout(3,2,'TileSpacing','tight','Padding','compact')

% colors
cols = [0.4660, 0.6740, 0.1880;0.4940, 0.1840, 0.5560;0.6350, 0.0780, 0.1840;...
    0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980];

%% Example patient, full spike rate
nexttile
rate = out(ex_p).rate;
cblock = out(ex_p).change_block;
% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(out(ex_p).unchanged_labels);
rate(ekg,:) = [];
run_dur = out(ex_p).run_dur;
block_dur = out(ex_p).block_dur;
rate = rate/run_dur;
rate = nanmean(rate,1); % sum across channels
times = (1:length(rate)) * block_dur;

% Get nan blocks
nan_blocks = find(isnan(nanmean(out(ex_p).rate,1)));

% plot
plot(times,rate,'k','linewidth',1)
ylabel('Spikes/elec/min')
 xlabel('Hour')
hold on

% Get surround times, starting with first non nan
surround = all_surrounds(is);
[pre,post] = get_surround_times(out(ex_p).rate,cblock,surround);
pre = pre*block_dur;
post = post*block_dur;
xlim([0 length(rate)*block_dur]);

set(gca,'fontsize',15)
yl = ylim;
new_yl = [yl(1) 1.3*(yl(2)-yl(1))];
top = yl(1) + 1.1*(yl(2)-yl(1));
ybar = yl(1) + 1.11*(yl(2)-yl(1));
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
plot([pre(1) pre(end)],[ybar ybar],'color',cols(1,:),'linewidth',4);
plot([post(1) post(end)],[ybar ybar],'color',cols(2,:),'linewidth',4);

xl = xlim;
cblock_fig_units = axescoord2figurecoord(cblock,nan);
ytext_fig_units = axescoord2figurecoord(ytext,nan);
annotation('textarrow',[0.2 cblock_fig_units-0.02],...
    [0.935 0.86],'String','Revision','color',cols(3,:),...
    'fontsize',15,'linewidth',2);

%{
text(pre(1),ytext,'Pre','color',cols(1,:),'fontsize',15,...
    'HorizontalAlignment','right')
text(post(end),ytext,'Post','color',cols(2,:),'fontsize',15,...
    'HorizontalAlignment','left')
    %}

pre_text_fig_units = axescoord2figurecoord(pre(1),nan);
pre_bar_fig_units = axescoord2figurecoord((pre(1)+pre(end))/2,nan);
post_text_fig_units = axescoord2figurecoord(post(end),nan);
post_bar_fig_units = axescoord2figurecoord((post(1)+post(end))/2,nan);

annotation('textarrow',[pre_text_fig_units-0.04 pre_text_fig_units-0.01],...
    [0.94 0.93],'String','Pre','color',cols(1,:),...
    'fontsize',15,'linewidth',2);

annotation('textarrow',[post_bar_fig_units+0.02 post_bar_fig_units-0.01],...
    [0.94 0.93],'String','Post','color',cols(2,:),...
    'fontsize',15,'linewidth',2);

xl = xlim;
yl = ylim;
%text(xl(1),yl(2),sprintf('Patient %d',ex_p),'fontsize',15,'VerticalAlignment','Top')
%xticklabels([])

%% Spike rate change across electrodes
nexttile
curr_ov = squeeze(all_ov(im,is,:,:));
plot([1 2],curr_ov','linewidth',2,'color','k')
hold on
xlim([0 3])
xticks([1 2])
xticklabels({'Pre','Post'})
% get the current tick labeks
ticklabels = get(gca,'XTickLabel');
% prepend a color for each tick label
ticklabels_new = cell(size(ticklabels));
for i = 1:length(ticklabels)
    ticklabels_new{i} = sprintf('\\color[rgb]{%f, %f, %f}%s',...
        cols(i,1),cols(i,2),cols(i,3),ticklabels{i});
end
% set the tick labels
set(gca, 'XTickLabel', ticklabels_new);
ylabel('Spikes/elec/min')

% Add stats
yl = ylim;
new_yl = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
ybar = yl(1) + 1.07*(yl(2)-yl(1));
yp = yl(1) + 1.13*(yl(2)-yl(1));
ylim(new_yl);
plot([1 2],[ybar ybar],'k')
text(1.5,yp,sprintf('%s',get_asterisks(main_stats(im,is,1),1)),...
    'horizontalalignment','center','fontsize',15);
set(gca,'fontsize',15)

% Results sentence
curr_rate_pvals = squeeze(all_rate_pval(is,:,im));
if any(curr_rate_pvals < 0.05/length(curr_rate_pvals))
    fprintf('\nSurprise significant individual patient rate result\n');
else
    fprintf(['\nNo individual patient had a larger peri-revision change in'...
        ' overall spike rates than that observed at randomly chosen time periods'...
        ' (Monte Carlo test with Bonferroni correction).\n']);
end

% Results Sentence
fprintf(['\nThere was no consistent difference across patients in the'...
    ' pre- (M = %1.1f, SD = %1.1f spikes/min) and post-revision'...
    ' (M = %1.1f, SD = %1.1f spikes/min) spike rate (paired t-test,'...
    ' t(%d) = %1.1f, %s)\n'],...
    mean(curr_ov(:,1)), std(curr_ov(:,1)),mean(curr_ov(:,2)),std(curr_ov(:,2)),...
    main_stats(im,is,3),main_stats(im,is,2),get_p_text(main_stats(im,is,1)))

%% Spike rate correlated with number of added electrodes
%{
nexttile
curr_ov = squeeze(all_ov(im,is,:,:));
ov_change = (curr_ov(:,2)-curr_ov(:,1))./(curr_ov(:,1));
plot(all_added_elecs(:,1),ov_change,'ko','markersize',15,'linewidth',2)
xlabel('# Added electrodes');
ylabel(sprintf('Rel. spike rate \\Delta'))
set(gca,'fontsize',15)
ylim([-1 2.5])
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('r = %1.2f\n%s',...
    added_stats(im,is,1),get_p_text(added_stats(im,is,2))),...
    'horizontalalignment','left',...
    'verticalalignment','top','fontsize',15)

% Results sentence
fprintf(['\nThere was a significant positive correlation across patients '...
    'between the number of electrodes added and the relative change in '...
    'spike rate (r = %1.2f, %s).\n'],added_stats(im,is,1),...
    get_p_text(added_stats(im,is,2)));
%}

%% Spike stability and node strength stability for main surround
% Loop over metrics (spikes vs node strength)
for im = 1:n_metrics
    
    metric = all_metrics{im};
    
    nexttile
    
    true_r = squeeze(all_all_true_r(is,:,im));
    mc_r = squeeze(all_all_mc_r(is,:,:,im));
    p = all_all_p(is,im);

    plot(1+0.05*rand(npts,1),true_r,'o','markersize',15,'linewidth',2,...
        'color',cols(4,:))
    hold on
    plot(2+0.05*rand(npts,1),mean(mc_r,2),'o','markersize',15,'linewidth',2,...
        'color',cols(5,:))
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    plot([1 2],[yl(1)+0.83*(yl(2)-yl(1)) yl(1)+0.83*(yl(2)-yl(1))],'k','linewidth',2)
    text(1.5,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(p,1),...
        'horizontalalignment','center','fontsize',15)
    xlim([0.5 2.5])
    xticks([1 2])
    xticklabels({'True','Monte Carlo'})
    if im == 1
        ylabel({'Spike stability'})
    else
        ylabel({'NS stability'})
    end
    set(gca,'fontsize',15)

    if im == 1
        
        % Results sentence
        curr_pvals = squeeze(all_all_all_p(is,:,im));
        if any(curr_pvals < 0.05/length(curr_pvals))
            fprintf('\nSurprise significant individual patient rho result\n');
        else
            fprintf(['\nNo individual patient had a larger peri-revision change in'...
                ' spike rate distribution than that observed at randomly chosen time periods'...
                ' (Monte Carlo test with Bonferroni correction).\n']);
        end
        
        fprintf(['\nFor the primary surround period of 24 hours, the spike stability '...
            'aggregated across patients (M = %1.2f, SD = %1.2f) was no different '...
            'from chance (Monte Carlo with Fisher''s method: %s) (Figure 4C).\n'],...
            mean(true_r),std(true_r),pretty_p_text(p));
    elseif im == 2
        
        % Results sentence
        curr_pvals = squeeze(all_all_all_p(is,:,im));
        if any(curr_pvals < 0.05/length(curr_pvals))
            fprintf('\nSurprise significant individual patient rho result\n');
        else
            fprintf(['\nNo individual patient had a larger peri-revision change in'...
                ' node strength distribution than that observed at randomly chosen time periods'...
                ' (Monte Carlo test with Bonferroni correction).\n']);
        end
        
        fprintf(['\n The group node strength stability in the surround period of 24 hours '...
            '(M = %1.2f, SD = %1.2f) was also no different from '...
            'chance (Monte Carlo with Fisher''s method: %s) (Figure 4D).\n'],...
            mean(true_r),std(true_r),pretty_p_text(p));
    end
    
    
end

%% Spike/node strength stability for all surrounds
for im = 1:n_metrics
    
    nexttile
    
    th = errorbar((1:n_surrounds)'-0.2,mean(all_all_true_r(:,:,im),2),std(all_all_true_r(:,:,im),[],2),...
        'o','linewidth',2,'markersize',10,'color',cols(4,:));
    hold on
    mch = errorbar((1:n_surrounds)'+0.2,mean(all_all_mc_r(:,:,:,im),[2 3]),std(all_all_mc_r(:,:,:,im),[],[2 3]),...
        'o','linewidth',2,'markersize',10,'color',cols(5,:));
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    for is = 1:n_surrounds
        plot([is-0.2,is+0.2],[yl(1)+0.8*(yl(2)-yl(1)) yl(1)+0.8*(yl(2)-yl(1))],...
            'k','linewidth',2)
        text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_all_p(is,im),n_surrounds),...
        'horizontalalignment','center','fontsize',15)
    end
    
    xlabel('Peri-revision period (hours)')
    xticks(1:n_surrounds)
    xticklabels(all_surrounds)
    xlim([0 n_surrounds+1])
    if im == 1
        ylabel({'Spike stability'})
        legend([th,mch],{'True','Monte Carlo'},'location','southeast','fontsize',15)
    else
        ylabel({'NS stability'})
    end
    set(gca,'fontsize',15)
    
    % Text
    fprintf(['\nExamining other peri-revision durations, no '...
            'duration had a spike stability significantly '...
            'lower than chance.\n']);
    %{
    if im == 1
        fprintf(['\nExamining other peri-revision surround durations, only the 6-hour '...
            'peri-revision surround duration had a spike stability significantly '...
            'lower than chance (spike stability M = %1.2f, SD = %1.2f, Monte Carlo %s).\n'],...
            mean(all_all_true_r(1,:,im)), std(all_all_true_r(1,:,im)),...
            get_p_text(all_all_p(1,im)));
    elseif im == 2
        fprintf(['\nExamining other surround durations, two '...
            'peri-revision surround durations (specifically ']); 
        sig_idx = find(all_all_p(:,im) < 0.05);
        for isig = 1:length(sig_idx)
            which_sig = sig_idx(isig);
            fprintf('%d',all_surrounds(which_sig));
            if isig < length(sig_idx) -1
                fprintf(', ');
            elseif isig == length(sig_idx) - 1
                fprintf(', and ')
            else
                fprintf(' hours');
            end
        end
        fprintf([') had a node strength stability significantly '...
            'lower than chance (see Figure 3 and Supplemental Table 1 for statistics).\n']);
    end
    %}

    
end

%% Annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.48 0.9 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.58 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.48 0.58 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.26 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
annotation('textbox',[0.48 0.26 0.1 0.1],'String','F','fontsize',20,'linestyle','none')
%annotation('textbox',[0.49 0.18 0.1 0.1],'String','G','fontsize',20,'linestyle','none')

%% Also save these to a specific range in the Supplemental Table 1
% rate stats
agg_rate_T = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.1f, %s',x,y,pretty_p_text(z)),...
    squeeze(main_stats(1,:,3))',squeeze(main_stats(1,:,2))',...
    squeeze(main_stats(1,:,1))','UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));
%{
agg_corr_T = cell2table(arrayfun(@(x,y) sprintf('r = %1.2f, %s',x,pretty_p_text(y)),...
    squeeze(added_stats(1,:,1))',squeeze(added_stats(1,:,2))','UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));
    %}
    
% spike stability    
agg_spike_T = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),all_all_p(:,1),...
    'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));

% node strength stability
agg_ns_T = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),all_all_p(:,2),...
    'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));

writetable(agg_rate_T,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','B2:B12','WriteVariableNames',false)
%writetable(agg_corr_T,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','E2:E12','WriteVariableNames',false)

writetable(agg_spike_T,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','C2:C12','WriteVariableNames',false)
writetable(agg_ns_T,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','D2:D12','WriteVariableNames',false)

%% Save fig
fname = 'new_rate';
print(gcf,[main_spike_results,fname],'-depsc');

end