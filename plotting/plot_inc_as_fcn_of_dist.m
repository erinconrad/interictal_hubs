function plot_inc_as_fcn_of_dist(rel_increase,abs_increase,dist_spikey,...
    spikey_labels,name,results_folder,run_dur,big_inc_labels)

%% Parameters
do_cat = 0;
do_spear = 1;
do_abs = 1;
do_ttest = 1;

if do_abs
    rate = abs_increase;
else
    rate = rel_increase;
end

if do_cat
    %% Prep output folder
    out_folder = [results_folder,'distance_cat/'];
    if ~exist(out_folder,'dir')
        mkdir(out_folder);
    end
    
    %% Are the big increase electrodes closer to closest added electrodes than others?
    all_chs = 1:length(dist_spikey);
    [is_big_inc] = ismember(spikey_labels,big_inc_labels);
    dist_inc = dist_spikey(is_big_inc);
    dist_no_inc = dist_spikey(~is_big_inc);
    
    if 0
        table(spikey_labels(is_big_inc))     
    end
    
    figure
    set(gcf,'position',[440   468   604   330])
    xpos_inc = 0.05*rand(length(dist_inc),1)+1;
    xpos_no_inc = 0.05*rand(length(dist_no_inc),1)+2;
    plot(xpos_inc,dist_inc,'o','color',[1 1 1])
    hold on
    text(xpos_inc,dist_inc,spikey_labels(is_big_inc),'fontsize',20,...
        'horizontalalignment','center')
    plot(xpos_no_inc,dist_no_inc,'o','color',[1 1 1])
    text(xpos_no_inc,dist_no_inc,spikey_labels(~is_big_inc),'fontsize',20,...
        'horizontalalignment','center')
    xticks([1 2])
    xticklabels({'Big spike increase','No big increase'})
    ylabel({'Distance from nearest','added electrode (mm)'})
    set(gca,'fontsize',20)
    xlim([0.5 2.5])
    
    %% Two sample T-test?
    if sum(is_big_inc) == 0
        title(sprintf('%s',name))
    else
        if do_ttest
            [~,pval] = ttest2(dist_inc,dist_no_inc);
            title(sprintf('%s two-sample t-test p = %1.3f',name,pval))
        else
            pval = ranksum(dist_inc,dist_no_inc);
            title(sprintf('%s Wilcoxon rank sum p = %1.3f',name,pval))
        end
    end
    
    %% save
    print(gcf,[out_folder,name,'_distance_cat'],'-dpng');
    close(gcf)
    
else

    %% Prep output folder
    out_folder = [results_folder,'distance/'];
    if ~exist(out_folder,'dir')
        mkdir(out_folder);
    end


    %% correlation between rate of increase and distance to closest added electrode
    if do_abs
        rate = rate/run_dur;
    end

    figure
    set(gcf,'position',[440 451 1001 347])
    plot(rate,dist_spikey,'o','color',[1 1 1]);
    hold on;
    text(rate,dist_spikey,spikey_labels,'horizontalalignment','center','fontsize',20)

    if do_abs
        xlabel('Change in spike rate after revision (spikes/min)')
    else
        xlabel('Relative increase in spike rate after revision')
    end
    ylabel({'Distance from nearest','added electrode (mm)'});
    if do_spear
        try
            [r,pval] = corr(rate,dist_spikey,'Type','Spearman','rows','complete');
        catch
            r = nan; pval = nan;
        end
        title(sprintf('%s Spearman correlation rho = %1.2f p = %1.3f',name,r,pval))
    else
        [r,pval] = corr(rate,dist_spikey,'Type','Pearson');
        title(sprintf('%s Pearson correlation r = %1.2f p = %1.3f',name,r,pval))
    end
    set(gca,'fontsize',20)

    print(gcf,[out_folder,name,'_distance'],'-dpng');
    close(gcf)

end



end