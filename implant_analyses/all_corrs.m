function all_corrs(whichPts,saved_out)

%% Parameters
surround = 24*0.5;
nb = 1e4;

which_resps = {'rel_rate','ns_rel'};
which_preds = {'dist','ns','cosi'};

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


figure
set(gcf,'position',[50 547 1000 800])
tiledlayout(2,3,'TileSpacing','compact','padding','compact')

% Loop over responses
for r = 1:length(which_resps)
    
    which_resp = which_resps{r};
    
    % Loop over predictors
    for p = 1:length(which_preds)
        
        which_pred = which_preds{p};
        
        % Initialize arrays of summary stats
        all_zs = nan(length(whichPts),7);
        all_mc_z = nan(length(whichPts),nb);
        
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
            
            % this is to make sure I am getting a reasonable z-score from my fisher
            % transformation.
            if abs(pval2-alt_pval2) > 0.05
                error('oh nos');
            end
            
            % Fill up all z's with info
            all_zs(i,:) = [z,z_score,rho,pval,n,pval2,alt_pval2];
            
        end
        
        % Stouffer's z-score
        z = nansum(all_zs(:,2))/sqrt(sum(~isnan(all_zs(:,2))));
        pval = 2*normcdf(-abs(z));
        r = tanh(nansum(all_zs(:,1).*all_zs(:,5))./nansum(all_zs(:,5)));
        
        % Get r back for MC Zs
        mc_r = nan(nb,1);
        for b = 1:nb
            mc_r(b) = tanh(nansum(all_mc_z(:,b).*all_zs(:,5))./nansum(all_zs(:,5)));
        end
        
        % Count number of mc_r's as or more extreme than true r
        num_more_sig = sum(abs(mc_r)>=abs(r));
        p_mc_agg = (num_more_sig + 1)/(nb+1);
        
        % Fisher's test to combine pvals
        % This tests whether, in aggregate, the revision time has a
        % disproportionately higher correlation relative to other times (regardless
        % of whether the direction is different across patients)
        X_2 = -2 * nansum(log(all_zs(:,4)));
        sum_p = 1-chi2cdf(X_2,2*sum(~isnan(all_zs(:,4))));
        
        % Plot aggregate statistics
        nexttile
        plot(all_zs(:,3),'o','markersize',20,'linewidth',2)
        hold on
        ylim([-1 1])
        plot(xlim,[0 0],'k--')
        yl = ylim;
        xl = xlim;
        text(xl(2),yl(2),...
            sprintf('Combined r = %1.2f\nMC p = %1.3f',...
            r,p_mc_agg),'fontsize',15,...
            'horizontalalignment','right',...
                'verticalalignment','top')
        xticklabels([])
        xlabel('Patient')
        title(sprintf('%s-%s',rtext,ptext))
        ylabel('Spearman correlation')
        set(gca,'fontsize',15)

        
    end  
end

print(gcf,[main_spike_results,sprintf('all_corrs_surround_%d',surround)],'-dpng')


end