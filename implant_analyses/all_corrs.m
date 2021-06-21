function all_corrs(whichPts,saved_out)

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
nb = 1e4;

n_surrounds = length(all_surrounds);
which_resps = {'rel_rate','ns_rel'};
%which_preds = {'dist','ns','cosi'};
which_preds = {'dist'};

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
all_p_simp = nan(length(all_surrounds),length(which_resps),length(which_preds),3);
all_r = nan(length(all_surrounds),length(which_resps),length(which_preds));
all_all_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));
all_all_mc_r = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));
all_all_p = nan(length(all_surrounds),length(which_resps),length(which_preds),length(whichPts));

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

            % Loop over patients
            for i = 1:length(whichPts) 

                rate = out(i).rate./out(i).run_dur;
                chLabels = out(i).unchanged_labels;
                cblock = out(i).change_block;
                ns = out(i).metrics.ns_norm;

                % Identify pre and post times
                [pre,post] = get_surround_times(rate,cblock,surround);

                ekg = identify_ekg_scalp(out(i).unchanged_labels);
                rate(ekg,:) = [];
                ns(ekg,:) = [];
                chLabels(ekg) = [];

                % Define predictor
                switch which_pred
                    case 'dist'
                        predictor = out(i).dist;
                        ptext = 'distance';  
                    case 'ns'
                        predictor = out(i).metrics.added_pc;
                        ptext = 'FC';  
                    case 'cosi'    
                        predictor = out(i).cosi;
                        ptext = 'CSI';
                end
                predictor(ekg) = [];

                % Define response
                switch which_resp
                    case 'rel_rate'
                        resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
                        rtext = 'Rate change';
                    case 'ns_rel'
                        resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./nanmean(ns(:,pre),2);
                        rtext = 'NS change';
                end


                % Monte carlo test
                [rho,pval,mc_rho] = mc_corr(rate,ns,predictor,...
                cblock,surround,nb,which_resp,'Spearman',0);

                % Fisher's R to z transformation on the original rhos
                n = sum(~isnan(resp) & ~isnan(predictor));
                [z,z_score,pval2] = fisher_transform(rho,n);
                [~,alt_pval2] = corr(resp,predictor,'Type','Spearman','rows','pairwise');

                % Also get the fisher transformed r-to-z for each mc rho
                mc_z = fisher_transform(mc_rho,n);
                all_mc_z(i,:) = mc_z;
                all_mc_r(i,:) = mc_rho;

                % this is to make sure I am getting a reasonable z-score from my fisher
                % transformation.
                if abs(pval2-alt_pval2) > 0.05
                    error('oh nos');
                end

                % Fill up all z's with info
                all_zs(i,:) = [z,z_score,rho,pval,n,pval2,alt_pval2];
                
                all_all_p(s,r,p,i) = pval;

            end

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
            all_p_mc(s,r,p) = p_mc_agg;
            all_p_simp(s,r,p,:) = [simp_p tstat df];
            all_r(s,r,p) = rho;
            all_all_r(s,r,p,:) = all_zs(:,3); % individual pt rhos
            all_all_mc_r(s,r,p,:) = mean(all_mc_r,2); % individual pt mean mc rhos
            
        end
    end
end

figure
set(gcf,'position',[50 547 700 500])
tt = tiledlayout(2,2,'TileSpacing','compact','padding','compact');


p = 1;

for r = 1:length(which_resps)

    if r == 1
        rtext = 'Rate change';
    else
        rtext = 'NS change';
    end

    
        
    if p == 1
        ptext = 'distance';
    elseif p == 2
        ptext = 'FC';  
    else
        ptext = 'CSI';  
    end

    % Plot aggregate statistics
    s = 1;
    nexttile
    plot(squeeze(all_all_r(s,r,p,:)),'o','markersize',15,'linewidth',2)
    hold on
    ylim([-1 1])
    xlim([1 length(whichPts)])
    plot(xlim,[0 0],'k--','linewidth',2)
    yl = ylim;
    xl = xlim;
    text(xl(2),yl(2),...
        sprintf('p = %1.3f',...
        all_p_mc(s,r,p)),'fontsize',15,...
        'horizontalalignment','right',...
            'verticalalignment','top')
    %{
    text(xl(2),yl(2),...
        sprintf('Combined r = %1.2f, p = %1.3f\nMC p = %1.3f',...
        all_r(s,r,p),all_p_simp(s,r,p,1),all_p_mc(s,r,p)),'fontsize',15,...
        'horizontalalignment','right',...
            'verticalalignment','top')
    %}
    xticklabels([])
    xlabel('Patient')
    ylabel(sprintf('%s-%s\ncorrelation',rtext,ptext))
    set(gca,'fontsize',15)
    
    % Plot for different surrounds
    nexttile
    
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
        text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_p_mc(is,r,p),1),...
        'horizontalalignment','center','fontsize',15)
    end
    xticks(1:n_surrounds)
    xticklabels(all_surrounds)
    xlabel('Peri-revision period (hours)')
    set(gca,'fontsize',15)
    
end

%% Annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.51 0.91 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.44 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.51 0.44 0.1 0.1],'String','D','fontsize',20,'linestyle','none')


print(gcf,[main_spike_results,sprintf('all_corrs_surround_%d',all_surrounds(1))],'-dpng')

%% Save table of p-values
names = ['all';names];
all_p_mc_spike = all_p_mc(:,1)';
all_p_mc_ns = all_p_mc(:,2)';

all_all_spike = [all_p_mc_spike;squeeze(all_all_p(:,1,1,:))'];
all_all_ns = [all_p_mc_ns;squeeze(all_all_p(:,2,1,:))'];

Tspike = cell2table(arrayfun(@(x) sprintf('%1.3f',x),all_all_spike','UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);

Tns = cell2table(arrayfun(@(x) sprintf('%1.3f',x),all_all_ns','UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);

writetable(Tspike,[main_spike_results,'spike_dist.csv'],'WriteRowNames',true)  
writetable(Tns,[main_spike_results,'ns_dist.csv'],'WriteRowNames',true)  


end