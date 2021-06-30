function all_corrs(whichPts,saved_out)

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
%all_surrounds = 12*[0.5,1];
main_surround = 3;
main_pred = 1;
main_resp = 1;
nb = 1e4;
ex = 1;

% Do fisher transformation on data? I set this to zero because of the
% potentially non-applicable assumptions of the transform. I get similar
% results either way and I don't think it is necessary given that my
% ultimate test for significance is a MC test. I keep it out for
% simplicity.
do_fisher = 0;

n_surrounds = length(all_surrounds);
which_resps = {'rel_rate','ns_rel'};
which_preds = {'dist','ns','cosi'};
%which_preds = {'dist'};

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

%% Main analysis
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

% Initialize stuff
all_p_mc = nan(length(all_surrounds),length(which_resps),length(which_preds));
all_p_mc_rho = nan(length(all_surrounds),length(which_resps),length(which_preds));

all_p_simp = nan(length(all_surrounds),length(which_resps),length(which_preds),3);
all_r = nan(length(all_surrounds),length(which_resps),length(which_preds));
all_all_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));
all_all_mc_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));
all_all_p = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));

%% Prep supplemental figure
figure
set(gcf,'position',[608 175 1371 514])
tiledlayout(2,5,'Padding','compact','tilespacing','compact')


% Loop over surrounds
for s = 1:length(all_surrounds)
    surround = all_surrounds(s);
    % Loop over responses
    for r = 1:length(which_resps)

        which_resp = which_resps{r};

        % Loop over predictors
        for p = 1:length(which_preds)

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
                ns = out(i).metrics.ns;

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
                        ptext = 'Distance';  
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
                end


                % Monte carlo test
                if p == main_pred
                    [rho,pval,mc_rho] = mc_corr(rate,ns,predictor,...
                    cblock,surround,nb,which_resp,'Spearman');
                else
                    rho = nan;
                    pval = nan;
                    mc_rho = nan;
                end
                

                % Fisher's R to z transformation on the original rhos
                n = sum(~isnan(resp) & ~isnan(predictor));
                [z,z_score,pval2] = fisher_transform(rho,n);
                [rho,alt_pval2] = corr(resp,predictor,'Type','Spearman','rows','pairwise');

                % Also get the fisher transformed r-to-z for each mc rho
                mc_z = fisher_transform(mc_rho,n);
                all_mc_z(i,:) = mc_z;
                all_mc_r(i,:) = mc_rho; % this is the MC rho I keep for main analysis
                all_true_r(i,:) = rho; % this is the true rho I keep for main analysis

                % this is to make sure I am getting a reasonable z-score from my fisher
                % transformation.
                if p == main_pred
                    if abs(pval2-alt_pval2) > 0.05
                        error('oh nos');
                    end
                end

                % Fill up all z's with info
                all_zs(i,:) = [z,z_score,rho,pval,n,pval2,alt_pval2];
                
                all_all_p(s,r,p,i) = pval; % these are the MC p values for individual patients
                
                %% Make supplemental figure
                if p == main_pred && s == main_surround && r == main_resp
                    
                    if alt_pval2 < 0.001
                        pvaltext = '< 0.001';
                    else
                        pvaltext = sprintf('= %1.3f',alt_pval2);
                    end
                    
                    nexttile
                    plot(predictor,resp,'o','linewidth',2)
                    ylabel(rtext)
                    xlabel(ptext)     
                    set(gca,'fontsize',15)
                    pause(0.2)
                    xl = xlim;
                    yl = ylim;
                    pause(0.2)
                    text(xl(2),yl(2),sprintf('Pt %d \\rho = %1.2f, p %s\nMC p = %1.3f',...
                        i,rho,pvaltext,pval),...
                        'HorizontalAlignment','Right','VerticalAlignment','Top',...
                        'fontsize',15)
                    
                end

            end

            %% Fisher transform combo
            % Get r back
            rho = tanh(nansum(all_zs(:,1).*all_zs(:,5))./nansum(all_zs(:,5)));
            [~,simp_p,~,stats] = ttest(all_zs(:,3));
            tstat = stats.tstat;
            df = stats.df;


            % Get r back for MC Zs
            mc_r = nan(nb,1);
            for b = 1:nb
                mc_r(b) = tanh(nansum(all_mc_z(:,b).*all_zs(:,5))./nansum(all_zs(:,5)));
            end

            % Count number of mc_r's as or more extreme than true r
            num_more_sig = sum(abs(mc_r)>=abs(rho));
            p_mc_agg = (num_more_sig + 1)/(nb+1);
            
            
            %% non fisher combo
            
            % Average true and MC rhos across patients (no fisher
            % transform)
            mean_true_rho = mean(all_true_r);
            mean_mc_rho = (mean(all_mc_r,1));
            
            % count how many MC rhos are as or more extreme
            num_more_sig_rho = sum(abs(mean_mc_rho)>=abs(mean_true_rho));
            p_mc_agg_rho = (num_more_sig_rho+1)/(nb+1);
            
            if 0
                figure
                plot(sort(mc_r),'o')
                hold on
                plot(xlim,[rho rho])
                if do_fisher == 1  
                    title(sprintf('p = %1.3f',p_mc_agg))
                else
                    title(sprintf('p = %1.3f',p_mc_agg_rho))
                end
                pause
                close(gcf)
            end

            % Fill up array
            all_p_mc(s,r,p) = p_mc_agg;
            all_p_mc_rho(s,r,p) = p_mc_agg_rho; % this is the one I keep for main analysi
            all_p_simp(s,r,p,:) = [simp_p tstat df];
            all_r(s,r,p) = rho;
            all_all_r(s,r,p,:) = all_zs(:,3); % individual pt rhos
            all_all_mc_r(s,r,p,:) = mean(all_mc_r,2); % individual pt mean mc rhos
            
        end
    end
end

% Save the supplemental figure
print(gcf,[main_spike_results,sprintf('supp_fig1_surround_%d',all_surrounds(main_surround))],'-dpng')


figure
%set(gcf,'position',[50 547 700 800])
set(gcf,'position',[50 329 687 468])
%{
tt = tiledlayout(3,2,'TileSpacing','compact','padding','compact');
tile_order = [1 3 5 2 4 6];
%}
tt = tiledlayout(2,2,'TileSpacing','compact','padding','compact');
tile_order = [1 3 2 4];

p = main_pred ; % distance
s = main_surround; 
count = 0;
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
    %{
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
    plot(1+0.05*rand(length(whichPts),1),squeeze(all_all_r(s,r,p,:)),'o','markersize',10,'linewidth',2)
    hold on
    plot(2+0.05*rand(length(whichPts),1),(squeeze(all_all_mc_r(s,r,p,:))),'o','markersize',10,'linewidth',2)
    ylim([-1 1])
    xlim([0.5 2.5])
    yl = ylim;
    ylim([yl(1) yl(1)+1.1*(yl(2)-yl(1))])
    yl = ylim;
    xl = xlim;
    plot([1 2],[yl(1)+0.7*(yl(2)-yl(1)) yl(1)+0.7*(yl(2)-yl(1))],...
            'k','linewidth',2)
    if do_fisher
        text(1.5,yl(1)+0.8*(yl(2)-yl(1)),get_asterisks(all_p_mc(1,r,p),1),...
        'horizontalalignment','center','fontsize',20)
        text(xl(2),yl(2),...
            sprintf('Combined r = %1.2f\nMC p = %1.3f',...
            all_r(s,r,p),all_p_mc(s,r,p)),'fontsize',15,...
            'horizontalalignment','right',...
                'verticalalignment','top')
    else
        text(1.5,yl(1)+0.8*(yl(2)-yl(1)),get_asterisks(all_p_mc_rho(1,r,p),1),...
        'horizontalalignment','center','fontsize',20)
    
        text(xl(2),yl(2),...
        sprintf('Combined r = %1.2f\nMC p = %1.3f',...
        mean(all_all_r(s,r,p,:),4),all_p_mc_rho(s,r,p)),'fontsize',15,...
        'horizontalalignment','right',...
            'verticalalignment','top')
    end
    
    %
    
    %}
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
    mch = errorbar((1:n_surrounds)'+0.2,mean(all_all_mc_r(:,r,p,:),4),std(all_all_mc_r(:,r,p,:),[],4),...
        'o','linewidth',2,'markersize',10);
        
    
    
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    for is = 1:n_surrounds
        plot([is-0.2,is+0.2],[yl(1)+0.8*(yl(2)-yl(1)) yl(1)+0.8*(yl(2)-yl(1))],...
            'k','linewidth',2)
        if do_fisher
            text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_p_mc(is,r,p),1),...
            'horizontalalignment','center','fontsize',15)
        else
            text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_p_mc_rho(is,r,p),1),...
            'horizontalalignment','center','fontsize',15)
        end
    end
    xticks(1:n_surrounds)
    xticklabels(all_surrounds)
    xlabel('Peri-revision period (hours)')
    ylabel(sprintf('%s-%s\ncorrelation',rtext,ptext))
    set(gca,'fontsize',15)
    
end

%% Annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.51 0.9 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.46 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.51 0.46 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
%annotation('textbox',[0 0.58 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
%annotation('textbox',[0.51 0.58 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
%annotation('textbox',[0 0.24 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
%annotation('textbox',[0.51 0.24 0.1 0.1],'String','F','fontsize',20,'linestyle','none')


print(gcf,[main_spike_results,sprintf('all_corrs_surround_%d',all_surrounds(main_surround))],'-dpng')

%% Save table of p-values
%{
p = 1;
names = ['all';names];
if do_fisher
    all_p_mc_spike = all_p_mc(:,1)';
    all_p_mc_ns = all_p_mc(:,2)';
else
    all_p_mc_spike = all_p_mc_rho(:,1,1)';
    all_p_mc_ns = all_p_mc_rho(:,2,1)';
end

all_all_spike = [all_p_mc_spike;squeeze(all_all_p(:,1,1,:))'];
all_all_ns = [all_p_mc_ns;squeeze(all_all_p(:,2,1,:))'];

Tspike = cell2table(arrayfun(@(x) sprintf('%1.3f',x),all_all_spike,'UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);

Tns = cell2table(arrayfun(@(x) sprintf('%1.3f',x),all_all_ns,'UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);

dtext = which_preds{1};
writetable(Tspike,[main_spike_results,'spike_',dtext,'.csv'],'WriteRowNames',true)  
writetable(Tns,[main_spike_results,'ns_',dtext,'.csv'],'WriteRowNames',true)  
%}

%% Also table of simple p-values
spike_p_simp = squeeze(all_p_simp(:,1,:,1));
ns_p_simp = squeeze(all_p_simp(:,2,:,1));

Tspike_simp = cell2table(arrayfun(@(x) sprintf('%1.3f',x),spike_p_simp,'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'VariableNames',which_preds);

Tns_simp = cell2table(arrayfun(@(x) sprintf('%1.3f',x),ns_p_simp,'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'VariableNames',which_preds);

writetable(Tspike_simp,[main_spike_results,'spike_simp_','.csv'],'WriteRowNames',true)  
writetable(Tns_simp,[main_spike_results,'ns_simp_','.csv'],'WriteRowNames',true)  

%% Sentences
fprintf(['\n\nOn a group level, the average correlation between relative spike rate change '...
    'and distance from the revision site was %1.2f, which was significantly '...
    'less than zero (t(%d) = %1.1f, p = %1.3f). Furthermore, this correlation '...
    'was stronger than expected for randomly chosen pseudo-revision times '...
    '(Monte Carlo p = %1.3f)\n'],...
    mean(all_all_r(main_surround,1,main_pred,:)),all_p_simp(main_surround,1,main_pred,3),...
    all_p_simp(main_surround,1,main_pred,2),all_p_simp(main_surround,1,main_pred,1),...
    all_p_mc_rho(main_surround,1,main_pred))

fprintf(['\n\nOn a group level, the average correlation between relative node strength change '...
    'and distance from the revision site was %1.2f, which was not significant '...
    't(%d) = %1.1f, p = %1.3f). This correlation '...
    'was not stronger than expected for randomly chosen pseudo-revision times '...
    '(Monte Carlo p = %1.3f)\n'],...
    mean(all_all_r(main_surround,2,main_pred,:)),all_p_simp(main_surround,2,main_pred,3),...
    all_p_simp(main_surround,2,main_pred,2),all_p_simp(main_surround,2,main_pred,1),...
    all_p_mc_rho(main_surround,2,main_pred))

end