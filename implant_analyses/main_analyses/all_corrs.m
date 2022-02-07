function all_corrs(whichPts,saved_out,out)

%{
This analysis tests whether electrodes closer to the revision site
experience a larger change in electrographic features peri-revision.

%}

%% Parameters
do_alt = 0;
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
main_surround = 3; %*******24 hours
main_pred = 1;
main_resp = 1;
nb = 1e1;  % change
do_save = 1;
do_buffer = 1;
type = 'Spearman';

%% Do combined probability test rather than more complicated method? (Should be 1)
do_comb_p = 1; % should be 1
do_fisher = 1; % Do fisher transformation on data? Should be 1
weighted_avg = 1; % weight by n-3 (doesn't matter, I no longer do this)

%% For supplemental table, do simple or mc p value? (should be 'simple')
which_p = 'simple'; % other is mc

n_surrounds = length(all_surrounds);
which_resps = {'rel_rate','ns_rel'};
which_preds = {'dist','ns','cosi'};
resp_text = {'relative spike rate','relative node strength'};

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
if do_alt
    main_spike_results = [results_folder,'main_spikes/alt/'];
else
    main_spike_results = [results_folder,'main_spikes/'];
end
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));


if isempty(whichPts)
    if do_alt
        whichPts = [20 103 105 106 107 108 35 109 110 94 97];
    else
        whichPts = [20 103 106 107 35 109 110 111 94 97];
    end
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
        if do_alt
            out(i) = alt_get_gdf_details(p);
        else
            out(i) = get_gdf_details(p);
        end
    end
    
    save([main_spike_results,'out'],'out');
    
end

%% Main analysis
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

%% Initialize stuff
% Group MC p-value (Fisher transformed or not)
all_p_mc = nan(length(all_surrounds),length(which_resps),length(which_preds));

% Group simple group p-values and associated stats
all_p_simp = nan(length(all_surrounds),length(which_resps),length(which_preds),3);

% Group rhos (either fisher tranformed or not)
all_r = nan(length(all_surrounds),length(which_resps),length(which_preds));

% Group MC rhos (either fisher transformed or not), 1 for each MC iteration
% and each patient
all_all_mc_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts),nb);

% individual rhos (not fisher transformed), 1 for each patient
all_all_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));

% Individual pt MC p-values
all_all_p = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));

% Occurrences of infinite relative change
inf_rel_change = nan(length(all_surrounds),length(which_resps),length(whichPts));

% Number old elecs
number_orig_elecs = nan(length(whichPts),1);

% Individual patient rhos (not fisher transformed), 1 for dist to depth, 1
% for dist to subdural (1 is depth, 2 is subdural)
all_all_r_type = nan(length(all_surrounds),length(which_resps),length(whichPts),2);
all_all_p_type = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts),2);
all_all_dist_resp = cell(length(all_surrounds),length(which_resps),length(whichPts),2); % distance and response for each electrode

% Distance to closest newest elecs
all_dist = cell(length(whichPts),1);

%% Prep supplemental figure
figure
set(gcf,'position',[608 175 1371 514])
tiledlayout(2,5,'Padding','compact','tilespacing','compact')
tilecount = 0;

% Loop over surrounds
for s = 1:length(all_surrounds)
    surround = all_surrounds(s);
    % Loop over responses
    for r = 1:length(which_resps) % relative rate change, relative node strength change

        which_resp = which_resps{r};

        % Loop over predictors
        for p = 1:length(which_preds) % distance, FC, cospike index

            which_pred = which_preds{p};

            % Initialize arrays of summary stats
            all_zs = nan(length(whichPts),7);
            all_mc_z = nan(length(whichPts),nb);
            all_mc_r = nan(length(whichPts),nb);
            all_true_r = nan(length(whichPts),1);

            % Loop over patients
            for i = 1:length(whichPts) 
                rate = out(i).rate./out(i).run_dur;
                chLabels = out(i).unchanged_labels;
                cblock = out(i).change_block;
   
                if isempty(out(i).metrics)
                    ns = nan(size(rate));
                else
                    ns = out(i).metrics.ns;
                end
                
                
                % Get file gap
                buffer = file_gaps(out(i).name);

                % Identify pre and post times
                [pre,post] = get_surround_times(rate,cblock,surround);

                % Remove ekg and scalp
                ekg = identify_ekg_scalp(out(i).unchanged_labels);
                rate(ekg,:) = [];
                ns(ekg,:) = [];
                chLabels(ekg) = [];

                % Define predictor
                switch which_pred
                    case 'dist'
                        predictor = out(i).dist;
                        ptext = 'Distance (mm)';  
                    case 'ns'
                        predictor = out(i).metrics.added_pc;
                        ptext = 'FC';  
                    case 'cosi'    
                        predictor = out(i).cosi;
                        ptext = 'CSI';
                end
                predictor(ekg) = []; % remove ekg
                
                % Remove those with zero predictor (a few funny electrodes
                % in HUP132 with poor localizations)
                zero_dist = predictor == 0;
                predictor(zero_dist) = [];
                rate(zero_dist,:) = [];
                ns(zero_dist,:) = [];
                chLabels(zero_dist) = [];

                % Define response
                switch which_resp
                    case 'rel_rate'
                        resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./abs(nanmean(rate(:,pre),2));
                        rtext = 'Rate change';
                    case 'ns_rel'
                        resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./abs(nanmean(ns(:,pre),2));
                        rtext = 'NS change';
                    case 'abs_rate'
                        resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2));
                        rtext = 'Rate change';
                    case 'ns_abs'
                        resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2));
                        rtext = 'NS change';
                end
                
                %% For revision stage, also correlate distance from different electrode types and spike rate change
                dist_depth = out(i).dist_depth;
                dist_subdural = out(i).dist_subdural;
                dist_depth(ekg) = []; % remove ekg
                dist_subdural(ekg) = []; % remove ekg
                dist_depth(zero_dist) = [];
                dist_subdural(zero_dist) = [];

                % Monte carlo test (I only do this for the distance
                % analysis)
                if p == main_pred
                    
                    % Get occurrences of infinite relative change
                    inf_rel_change(s,r,i) = sum(resp == inf);
                    number_orig_elecs = length(resp);
                    
                    % Distance to closest new elec
                    all_dist{i} = predictor;

                    % MC test
                    [rho,pval,mc_rho] = mc_corr(rate,ns,predictor,...
                    cblock,surround,nb,which_resp,type,buffer,do_buffer);
                
                    %% For revision stage, also correlate distance from different electrode types and spike rate change
                    %{
                    [rho_depth,pval_depth,mc_rho_depth] = mc_corr(rate,ns,dist_depth,...
                        cblock,surround,nb,which_resp,type,buffer,do_buffer);
                    

                    [rho_subdural,pval_subdural,mc_rho_subdural] = mc_corr(rate,ns,dist_subdural,...
                        cblock,surround,nb,which_resp,type,buffer,do_buffer);
                    %}
                    rho_depth = corr(resp,dist_depth,'Type',type,'rows','pairwise');
                    rho_subdural = corr(resp,dist_subdural,'Type',type,'rows','pairwise');
                    all_all_r_type(s,r,i,:) = [rho_depth,rho_subdural];
                    
                    % get the response and the distance from nearest depth
                    % and subdural
                    all_all_dist_resp{s,r,i,1} = [resp,dist_depth];
                    all_all_dist_resp{s,r,i,2} = [resp,dist_subdural];
                    
                    %all_all_mc_r_type(s,r,p,i,:) = [mc_rho_depth,mc_rho_subdural];
                    %all_all_p_type(s,r,p,i,:) = [pval_depth,pval_subdural];
                    
                else
                    rho = corr(resp,predictor,'Type',type,'rows','pairwise');
                    pval = nan;
                    mc_rho = nan;
                end
                
                
                

                % Fisher's R to z transformation on the original rhos
                n = length(chLabels); 
                % this n isn't exactly the number involved in the
                % correlation because some electrodes have nans in some
                % blocks, but this number changes throughout the time
                % periods and so this is a fair approximation for the
                % purpose of weighting across patients
                
                [z,~,~] = fisher_transform(rho,n);
                [rho,alt_pval2] = corr(resp,predictor,'Type',type,'rows','pairwise');

                % Also get the fisher transformed r-to-z for each mc rho
                mc_z = fisher_transform(mc_rho,n);
                all_mc_z(i,:) = mc_z; % MC fisher r-to-z transformation (which I will compare against true z)
                all_mc_r(i,:) = mc_rho; % MC rho 
                all_true_r(i) = rho; % True rho

                % Fill up all z's with info
                all_zs(i,:) = [z,nan,rho,pval,n,nan,alt_pval2];
                
                all_all_p(s,r,p,i) = pval; % these are the MC p values for individual patients
                
                %% Make supplemental figure for primary conditions
                % This shows individual patient correlations
                if p == main_pred && s == main_surround && r == main_resp
                    
                    % alt_pval2 is the non MC p value for the individual
                    % patient correlation
                    pvaltext = pretty_p_text(alt_pval2);
                    mcpvaltext = pretty_p_text(pval);
                    
                    nexttile
                    tilecount = tilecount + 1;
                    plot(predictor,resp,'o','linewidth',2)
                    hold on
                                        
                    % add infinite values at some highest value
                    inf_resp = resp == inf;
                    if sum(inf_resp) > 0
                        yl = ylim;
                        old_yl_2 = yl(2);
                        new_yl_2 = yl(1) + 1.1*(yl(2)-yl(1));
                        midpoint = mean([yl(2) new_yl_2]);
                        new_yl_3 = yl(1) + 1.13*(yl(2)-yl(1));
                        col = [0.8500, 0.3250, 0.0980];
                        plot(predictor(inf_resp),new_yl_2,'*','color',col,'linewidth',2)
                        ylim([yl(1) new_yl_3])
                    else
                        yl = ylim;
                        old_yl_2 = yl(2);
                    end
                    if tilecount == 1 || tilecount == 6
                        ylabel('Relative rate change')
                    end
                    if tilecount > 5
                        xlabel(ptext)     
                    end
                    set(gca,'fontsize',15)
                    pause(0.2) % delete at your own risk
                    xl = xlim;
                    %yl = ylim;
                    pause(0.2)
                    text(xl(2),old_yl_2,sprintf('Pt %d\n\\rho = %1.2f, %s\nMC %s',...
                        i,rho,pvaltext,mcpvaltext),...
                        'HorizontalAlignment','Right','VerticalAlignment','Top',...
                        'fontsize',15)
                    if sum(inf_resp) > 0
                        % change yaxis to clarify break
                        yticks([yl(1) old_yl_2 new_yl_2])
                        yticklabels({sprintf('%1.1f',yl(1)),sprintf('%1.1f',old_yl_2),...
                            '\infty'});
                        plot(xlim,[midpoint midpoint],'k--');
                    end
                end

            end

            
            
            %% Get simple (non-MC) group stats
            if do_fisher
                % One sample t-test on the z's
                [~,simp_p,~,stats] = ttest(all_zs(:,1));
            else
                % One sample t-test on the rho's
                [~,simp_p,~,stats] = ttest(all_true_r);
            end
            tstat = stats.tstat;
            df = stats.df;
            
            
            %% Get MC stats
            if do_comb_p
                
                % Fisher test to combine pvalues for ros
                if sum(isnan(all_all_p(s,r,p,:))) > 0 && p == 1, error('why'); end
                X_2 = -2 * sum(log(all_all_p(s,r,p,:)));
                sum_p = 1-chi2cdf(X_2,2*length(all_all_p(s,r,p,:)));

                p_mc_agg = sum_p;
            %{    
            elseif do_fisher
                % Fisher transform combo analysis
                % Get r back
                if weighted_avg
                    % this is tanh(sum(z*(n-3))/sum(n-3))
                    rho = tanh(nansum(all_zs(:,1).*(all_zs(:,5)-3))./nansum((all_zs(:,5)-3)));
                else
                    rho = tanh(nanmean(all_zs(:,1)));
                end


                % Get r back for MC Zs
                mc_r = nan(nb,1);
                for b = 1:nb
                    if weighted_avg
                        % same n's
                        mc_r(b) = tanh(nansum(all_mc_z(:,b).*(all_zs(:,5)-3))./nansum((all_zs(:,5)-3)));
                    else
                        mc_r(b) = tanh(nanmean(all_mc_z(:,b)));
                    end
                end

                % Count number of mc_r's as or more extreme than true r -
                % group MC analysis
                num_more_sig = sum(abs(mc_r)>=abs(rho));
                p_mc_agg = (num_more_sig + 1)/(nb+1);
                                
            else
                %Non fisher transform combo analysis
                
                % Average true and MC rhos across patients (no fisher
                % transform)
                rho = mean(all_zs(:,3));
                mc_r = (mean(all_mc_r,1));
                
                % Count number of mc_r's as or more extreme than true r -
                % group MC analysis
                num_more_sig = sum(abs(mc_r)>=abs(rho));
                p_mc_agg = (num_more_sig + 1)/(nb+1);
                %}
                
            end
            
                 
            
            if 0
                figure
                plot(sort(mc_r),'o')
                hold on
                plot(xlim,[rho rho])
                title(sprintf('p = %1.3f',p_mc_agg))
                pause
                close(gcf)
            end

            % Fill up array
            all_p_mc(s,r,p) = p_mc_agg; % mc analysis
            all_p_simp(s,r,p,:) = [simp_p tstat df]; % simple analysis
            all_r(s,r,p) = rho; % group rho (fisher or not)
            all_all_r(s,r,p,:) = all_zs(:,3); % individual pt rhos (not transformed)
            all_all_mc_r(s,r,p,:,:) = all_mc_r; % individual mc rhos, all iterations
            
        end
    end
end

% Save the supplemental figure
if do_save
    print(gcf,[main_spike_results,sprintf('supp_fig1_surround_%d',all_surrounds(main_surround))],'-depsc')
end

%% Say how many occurrences of inifinte relative change there are
s = main_surround;

    
fprintf(['\nThis occurred for an average of %1.1f (%1.1f%%) electrodes '...
    'across patients.\n'],...
    mean(inf_rel_change(s,1,:)),mean(inf_rel_change(s,1,:)./number_orig_elecs*100));
    

%% Say distance to closest new electrode
all_all_dist = [];
for i = 1:length(all_dist)
   all_all_dist = [all_all_dist;all_dist{i}]; 
end
min_dist = cellfun(@min,all_dist);
fprintf(['\nAcross all patients and all electrodes, the mean (SD) distance '...
    'to the revision site was %1.1f (%1.1f) mm. '...
    'Across patients, the original electrode closest to the revision site was '...
    'on average %1.1f (SD %1.1f) mm from the revision site\n'],...
    mean(all_all_dist),std(all_all_dist),...
    mean(min_dist),std(min_dist));

%% Do anatomy analysis
[output_for_plot,loc_names] = anatomy_analyses(whichPts,saved_out,out);

%% initialize main figure
figure
%set(gcf,'position',[50 547 700 800])
set(gcf,'position',[50 329 687 700])
%{
tt = tiledlayout(3,2,'TileSpacing','compact','padding','compact');
tile_order = [1 3 5 2 4 6];
%}
tt = tiledlayout(3,2,'TileSpacing','compact','padding','compact');
tile_order = [1 3 2 4];

% set predictor and surround
p = main_pred ; % distance
s = main_surround; % 24 hours
count = 0;
% Loop over rate and node strength
for r = 1:length(which_resps)
    which_resp = which_resps{r};
    if r == 1
        rtext = 'Rate change';
    else
        rtext = 'NS change';
    end

    
        
    if p == 1
        which_pred = which_preds{p};
        ptext = 'distance';
    elseif p == 2
        ptext = 'FC';  
    else
        ptext = 'CSI';  
    end
    
    %% Plot single pt example
    %{
    labels = out(ex).unchanged_labels;
    cblock = out(ex).change_block;
    rate = out(ex).rate;
    ekg = identify_ekg_scalp(labels);
    ns = out(ex).metrics.ns;
    
    % remove ekg
    rate(ekg,:) = [];
    ns(ekg,:) = [];
      
    % Surround
    surround = all_surrounds(s);
    [pre,post] = get_surround_times(rate,cblock,surround);
    
    % Define response
    switch which_resp
        case 'rel_rate'
            resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
            rtext = 'Rate change';
        case 'ns_rel'
            resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./nanmean(ns(:,pre),2);
            rtext = 'NS change';
    end
    
    % Define predictor
    switch which_pred
        case 'dist'
            predictor = out(ex).dist;
            ptext = 'distance';  
        case 'ns'
            predictor = out(ex).metrics.added_pc;
            ptext = 'FC';  
        case 'cosi'    
            predictor = out(ex).cosi;
            ptext = 'CSI';
    end
    predictor(ekg) = [];
    
    % Plot the correlation
    
    count = count+1;
    nexttile(tile_order(count));
    plot(predictor,resp,'o','markersize',10,'linewidth',2)
    xlabel('Distance (mm)')
    ylabel(rtext)
    yl = ylim;
    xl = xlim;
    set(gca,'fontsize',15)
    %pause
    text(xl(2),yl(2),...
        sprintf('r = %1.2f\nMC p = %1.2f',...
        (all_all_r(s,r,p,ex)),all_all_p(s,r,p,ex)),'fontsize',15,...
        'horizontalalignment','right',...
            'verticalalignment','top')
    
    %}

    %% Plot aggregate statistics
    count = count+1;
    nexttile(tile_order(count));
    % Plot the individual patient rho's (not transformed)
    plot(1+0.08*randn(length(whichPts),1),squeeze(all_all_r(s,r,p,:)),'o','markersize',10,'linewidth',2)
    hold on
    % plot the mean individual patient MCs (not transformed), averaged over
    % iterations
    errorbar(2+0.08*randn(length(whichPts),1),...
        (squeeze(mean(all_all_mc_r(s,r,p,:,:),5))),...
        (squeeze(std(all_all_mc_r(s,r,p,:,:),[],5))),...
        'o','markersize',10,'linewidth',2)
    ylim([-1 1])
    xlim([0.5 2.5])
    yl = ylim;
    ylim([yl(1) yl(1)+1.1*(yl(2)-yl(1))])
    yl = ylim;
    xl = xlim;
    plot([1 2],[yl(1)+0.75*(yl(2)-yl(1)) yl(1)+0.75*(yl(2)-yl(1))],...
            'k','linewidth',2)
    if strcmp(which_p,'simple')
        text(1.5,yl(1)+0.82*(yl(2)-yl(1)),get_asterisks(all_p_simp(s,r,p),1),...
        'horizontalalignment','center','fontsize',20)
    
        text(xl(2),yl(2),...
            sprintf('Combined r = %1.2f\n%s',...
            all_r(s,r,p),pretty_p_text(all_p_simp(s,r,p))),'fontsize',15,...
            'horizontalalignment','right',...
                'verticalalignment','top')
        
    elseif strcmp(which_p,'mc')
        text(1.5,yl(1)+0.82*(yl(2)-yl(1)),get_asterisks(all_p_mc(s,r,p),1),...
        'horizontalalignment','center','fontsize',15)
    
        text(xl(2),yl(2),...
            sprintf('Combined r = %1.2f\nMC %s',...
            all_r(s,r,p),pretty_p_text(all_p_mc(s,r,p))),'fontsize',15,...
            'horizontalalignment','right',...
                'verticalalignment','top')
    end
    
    xticks([1 2])
    xticklabels({'True','Monte Carlo'})
    ylabel(sprintf('%s-%s\ncorrelation',rtext,ptext))
    set(gca,'fontsize',15)
    
    %% Plot for different surrounds
    count = count+1;
    nexttile(tile_order(count));
    
    th = errorbar((1:n_surrounds)'-0.2,mean(all_all_r(:,r,p,:),4),std(all_all_r(:,r,p,:),[],4),...
    'o','linewidth',2,'markersize',10);
    hold on
    mch = errorbar((1:n_surrounds)'+0.2,mean(all_all_mc_r(:,r,p,:,:),[4 5]),std(all_all_mc_r(:,r,p,:,:),[],[4 5]),...
        'o','linewidth',2,'markersize',10);
        
    
    
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    for is = 1:n_surrounds
        plot([is-0.2,is+0.2],[yl(1)+0.8*(yl(2)-yl(1)) yl(1)+0.8*(yl(2)-yl(1))],...
            'k','linewidth',2)
        if strcmp(which_p,'simple')
            text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_p_simp(is,r,p),n_surrounds),...
            'horizontalalignment','center','fontsize',15)
        elseif strcmp(which_p,'mc')
            text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_p_mc(is,r,p),n_surrounds),...
            'horizontalalignment','center','fontsize',15)
        end
        
    end
    xlim([0 n_surrounds+1])
    xticks(1:n_surrounds)
    xticklabels(all_surrounds)
    xtickangle(45)
    xlabel('Peri-revision period (hours)')
    ylabel(sprintf('%s-%s\ncorrelation',rtext,ptext))
    set(gca,'fontsize',15)
    
end

%% do anatomy plots
for i = 1:2 % spikes and node strength
    nexttile
    for j = 1:4
        data = squeeze(output_for_plot(:,j,i));
        errorbar(j,mean(data),std(data),'o','linewidth',2,'markersize',10);
        hold on
    end
    p = friedman(output_for_plot(:,:,i),1,'off');
    if p < 0.05, error('woah'); end
    
    yl = ylim;
    ybar = yl(1) + 1.02*(yl(2)-yl(1));
    ytext = yl(1) + 1.1*(yl(2)-yl(1));
    ylnew = [yl(1) yl(1) + 1.15*(yl(2)-yl(1))];
    ylim(ylnew)
    
    plot([1 4],[ybar ybar],'k-','linewidth',2)
    text(2.5,ytext,'ns','horizontalalignment','center','fontsize',15)
    
    xlim([0 5])
    xticks(1:4)
    xticklabels(loc_names)
    xtickangle(30)
    if i == 1
        ylabel({'Relative spike','rate change'});
    else
        ylabel({'Relative node','strength change'});
    end
    set(gca,'fontsize',15)
end


%% Annotations
annotation('textbox',[0.03 0.90 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.54 0.90 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0.03 0.59 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.54 0.59 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0.03 0.28 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
annotation('textbox',[0.54 0.28 0.1 0.1],'String','F','fontsize',20,'linestyle','none')
%annotation('textbox',[0 0.58 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
%annotation('textbox',[0.51 0.58 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
%annotation('textbox',[0 0.24 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
%annotation('textbox',[0.51 0.24 0.1 0.1],'String','F','fontsize',20,'linestyle','none')

if do_save 
    print(gcf,[main_spike_results,sprintf('all_corrs_surround_%d',all_surrounds(main_surround))],'-depsc')
end

%% Save table of p-values
spike_p_simp = squeeze(all_p_simp(:,1,1,:));
ns_p_simp = squeeze(all_p_simp(:,2,1,:));

Tspike_simp = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.1f, %s',x,y,pretty_p_text(z)),...
    spike_p_simp(:,3),spike_p_simp(:,2),spike_p_simp(:,1),'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));

Tns_simp = cell2table(arrayfun(@(x,y,z) sprintf('t(%d) = %1.1f, %s',x,y,pretty_p_text(z)),...
    ns_p_simp(:,3),ns_p_simp(:,2),ns_p_simp(:,1),'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));



    %writetable(Tspike_simp,[main_spike_results,'spike_simp_','.csv'],'WriteRowNames',true)  
    %writetable(Tns_simp,[main_spike_results,'ns_simp_','.csv'],'WriteRowNames',true)  

% Also table of MC p-values
spike_p_MC = squeeze(all_p_mc(:,1,main_pred));
ns_p_MC = squeeze(all_p_mc(:,2,main_pred));

Tspike_MC = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),spike_p_MC,'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));

Tns_MC = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),ns_p_MC,'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));


    %writetable(Tspike_MC,[main_spike_results,'spike_MC','.csv'],'WriteRowNames',true)  
    %writetable(Tns_MC,[main_spike_results,'ns_MC','.csv'],'WriteRowNames',true)  
    
if do_save
    if strcmp(which_p,'simple')
        %% Also add to main supplemental table
        writetable(Tspike_simp,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','E2:E12','WriteVariableNames',false)
        writetable(Tns_simp,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','F2:F12','WriteVariableNames',false)
    elseif strcmp(which_pc,'mc')
    %% Also add to main supplemental table
        writetable(Tspike_MC,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','E2:E12','WriteVariableNames',false)
        writetable(Tns_MC,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','F2:F12','WriteVariableNames',false)
    end
end

%% Sentences

sp_p_main = squeeze(all_all_p(3,1,1,:));
ns_p_main = squeeze(all_all_p(3,2,1,:));

if any(sp_p_main < 0.05/length(sp_p_main))
    fprintf('\nSurprise significant individual patient spike correlation result\n');
else
    fprintf(['\nNo individual patient had a larger peri-revision correlation between relative spike rate change and'...
        ' distance from the implant revision site than that observed at randomly chosen time periods'...
        ' (Monte Carlo test with Bonferroni correction).\n']);
end

if any(ns_p_main < 0.05/length(ns_p_main))
    fprintf('\nSurprise significant individual patient rho result\n');
else
    fprintf(['\nNo individual patient had a larger peri-revision correlation between relative node strength change and'...
        ' distance from the implant revision site than that observed at randomly chosen time periods'...
        ' (Monte Carlo test with Bonferroni correction).\n']);
end

fprintf(['\n\nThe average correlation across patients between relative spike rate change '...
    'and distance from the revision site was %1.2f, which was not significantly '...
    'different from zero (t(%d) = %1.1f, p = %1.3f). Furthermore, this correlation '...
    'was not stronger than expected for randomly chosen pseudo-revision times '...
    '(Monte Carlo p = %1.3f)\n'],...
    mean(all_all_r(main_surround,1,main_pred,:)),all_p_simp(main_surround,1,main_pred,3),...
    all_p_simp(main_surround,1,main_pred,2),all_p_simp(main_surround,1,main_pred,1),...
    all_p_mc(main_surround,1,main_pred))

fprintf(['\n\nThe average correlation between relative node strength change '...
    'and distance from the revision site was %1.2f, which was not significant '...
    't(%d) = %1.1f, p = %1.3f). This correlation '...
    'was not stronger than expected for randomly chosen pseudo-revision times '...
    '(Monte Carlo p = %1.3f)\n'],...
    mean(all_all_r(main_surround,2,main_pred,:)),all_p_simp(main_surround,2,main_pred,3),...
    all_p_simp(main_surround,2,main_pred,2),all_p_simp(main_surround,2,main_pred,1),...
    all_p_mc(main_surround,2,main_pred))


% all_all_r has dimensions n_surround x n_response x n_predictor x
% n_patients
% and so all_all_r(main_surround,1,2,:) is main surround, spike rate, FC,
% all patients
% all_all_r(main_surround,1,3,:) is main surround, spike rate, cosi, all
% patients
fprintf(['\n\nExamining other measures of proximity to the revision site,'...
    ' on a group level there was no significant correlation between relative spike rate'...
    ' change and either functional connectivity (average rho = %1.2f, t(%d) = %1.1f, p = %1.2f)'...
    ' or co-spike index (average rho = %1.2f, t(%d) = %1.1f, p = %1.2f) with the revision site.\n\n'],...
    mean(all_all_r(main_surround,1,2,:)),all_p_simp(main_surround,1,2,3),...
    all_p_simp(main_surround,1,2,2),all_p_simp(main_surround,1,2,1),...
     mean(all_all_r(main_surround,1,3,:)),all_p_simp(main_surround,1,3,3),...
    all_p_simp(main_surround,1,3,2),all_p_simp(main_surround,1,3,1));

% all_all_r(main_surround,2,2,:) is main surorund, node strength, FC, all
% patients
fprintf(['\n\nExamining other measures of proximity to the revision site,'...
    ' on a group level there was no significant correlation between relative node strength'...
    ' change and either functional connectivity (average rho = %1.2f, t(%d) = %1.1f, p = %1.2f)'...
    ' or co-spike index (average rho = %1.2f, t(%d) = %1.1f, p = %1.2f) with the revision site.\n\n'],...
    mean(all_all_r(main_surround,2,2,:)),all_p_simp(main_surround,2,2,3),...
    all_p_simp(main_surround,2,2,2),all_p_simp(main_surround,2,2,1),...
     mean(all_all_r(main_surround,2,3,:)),all_p_simp(main_surround,2,3,3),...
    all_p_simp(main_surround,2,3,2),all_p_simp(main_surround,2,3,1));


%% For reviews, test whether original electrodes very close to depths have a bigger spike rate change than those close to subdurals
dthreshs = [10 20 30];
type_table = cell(length(all_surrounds),2,length(dthreshs));
type_text = {'depth','subdural'};
out_text = cell(2,1);

% loop over threshold distances
for id = 1:length(dthreshs)
    dthresh = dthreshs(id);
    for is = 1:length(all_surrounds)
        for ir = 1:2 % spike rate and node strength
            for ipatient = 2
                
                % get response and dist for depth and subdural-proximate
                % electrodes
                depth = squeeze(all_all_dist_resp{is,ir,ipatient,1});
                subdural = squeeze(all_all_dist_resp{is,ir,ipatient,2});
                
                % find those electrodes very close to added depths and subdurals,
                % respectively
                close_depth = depth(:,2) < dthresh;
                close_subdural = subdural(:,2) < dthresh;
                
                % difference in relative feature change between those
                % close to depths and those close to subdurals?
                [pval,~,stats] = ranksum(depth(close_depth,1),subdural(close_subdural,1));
                W = stats.ranksum;
                nt = sum(~(isnan(depth(close_depth,1))));
                ne = sum(~(isnan(subdural(close_subdural,1))));
                U1 = W - nt*(nt+1)/2;
                U2 = nt*ne-U1;
                U = min([U1,U2]);
                
                if id == 2 && is == main_surround
                    fprintf(['\nThe median %s change for original contacts within'...
                        ' %d mm of the closest added %s contacts was (%1.2f), which was not'...
                        ' significantly different from that for contacts within %d mm of the'...
                        ' closest added %s contacts (%1.2f) (Mann-Whitney test: U'...
                        '(N_depth_proximate = %d, N_subdural_proximate = %d) = %1.1f, %s).\n'],...
                        resp_text{ir},dthresh,type_text{1},nanmedian(depth(close_depth,1)),...
                        dthresh,type_text{2},nanmedian(depth(close_subdural,1)),nt,ne,U,simple_p_text(pval));
                   % out_text{ir} = ttext;
                end
                table_text = sprintf('U(Nd = %d, Ns = %d) = %1.1f, %s',nt,ne,U,simple_p_text(pval));
                type_table{is,ir,id} = table_text;

                if 0
                    figure
                    plot(1+randn(sum(close_depth),1)*0.05,depth(close_depth,1),'o')
                    hold on
                    plot(2+randn(sum(close_subdural),1)*0.05,subdural(close_subdural,1),'o')
                    title(sprintf('p = %1.3f',pval))
                    pause
                    close gcf
                end
                

            end
        end
    end
end

new_type_table = cell2table([{'10 mm threshold spike rate change','10 mm threshold node strength change'};...
    type_table(:,:,1);{'',''};{'20 mm threshold spike rate change','20 mm threshold node strength change'};...
    type_table(:,:,2);{'',''};{'30 mm threshold spike rate change','30 mm threshold node strength change'};type_table(:,:,3)]);
writetable(new_type_table,[main_spike_results,'Supplemental Table 4.xlsx'],'Range','A2:B39','WriteVariableNames',false)

%out_text{1}
%out_text{2}
%{
%figure
%tiledlayout(2,2)
type_table = nan(2,2,2); % spikes vs ns, patient, depth vs subdural
for ir = 1:2 % spikes and ns
    curr_r_type = squeeze(all_all_r_type(is,ir,ip,:,:));
    curr_p_type = squeeze(all_all_p_type(is,ir,ip,:,:));
    
    % find patients without any nans
    no_nans = ~any(isnan(curr_r_type),2);
    
    curr_r_type = curr_r_type(no_nans);
    curr_p_type = curr_p_type(no_nans);
    npts = size(curr_r_type,1);
    for p = 1:npts
        curr_curr = curr_r_type(p,:);
    end
end
%}
    
end