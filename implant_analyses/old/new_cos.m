function new_cos(post,cosi,abs_increase,spikey_idx,labels,name,results_folder,dist,rel_increase,do_rel)

%% Just do spikey
post = post(spikey_idx);
cosi = cosi(spikey_idx);
dist = dist(spikey_idx);
labels = labels(spikey_idx);
abs_increase = abs_increase(spikey_idx);
rel_increase = rel_increase(spikey_idx);

%% get co-spike index
%cosi = cos./post;
if sum(spikey_idx) == 0
    fprintf('\nNo spikey electrodes for %s, skipping...\n',name);
    return
end

figure
set(gcf,'position',[60 451 1350 355])

%% Are electrodes close to the newly added electrodes more likely to co-spike with them?
subplot(1,3,1)
plot(dist,cosi,'o','color',[1 1 1])
hold on
text(dist,cosi,labels,'horizontalalignment','center','fontsize',20)
xlabel('Distance from new electrodes')
ylabel('Co-spike index')
if isempty(cosi)
    rho = nan; pval = nan;
else
    [rho,pval] = corr(dist,cosi,'Type','Spearman','rows','pairwise');
end
title(sprintf('%s rho = %1.2f p = %1.3f',name,rho,pval))
set(gca,'fontsize',20);

%% Are electrodes close to newly added electrodes more likely to have a change in spike rate
subplot(1,3,2)
if do_rel
    plot(dist,rel_increase,'o','color',[1 1 1])
    hold on
    text(dist,rel_increase,labels,'horizontalalignment','center','fontsize',20)
    ylabel({'Relative change in spike rate','after revision'})
else
    plot(dist,abs_increase,'o','color',[1 1 1])
    hold on    
    text(dist,abs_increase,labels,'horizontalalignment','center','fontsize',20)
    ylabel({'Change in spike rate','after revision (spikes/min)'})
end
xlabel('Distance from new electrodes')
if do_rel
    [rho,pval] = corr(rel_increase,dist,'Type','Spearman','rows','pairwise');
else
    [rho,pval] = corr(abs_increase,dist,'Type','Spearman','rows','pairwise');
end

title(sprintf('%s rho = %1.2f p = %1.3f',name,rho,pval))
set(gca,'fontsize',20);

%% Are electrodes that co-spike with the new electrodes more likely to have a rate change
subplot(1,3,3)
if do_rel
    plot(cosi,rel_increase,'o','color',[1 1 1])
    hold on
    text(cosi,rel_increase,labels,'horizontalalignment','center','fontsize',20)
    ylabel({'Relative change in spike rate','after revision'})
else
    plot(cosi,abs_increase,'o','color',[1 1 1])
    hold on
    text(cosi,abs_increase,labels,'horizontalalignment','center','fontsize',20)
    ylabel({'Change in spike rate','after revision (spikes/min)'})
end
xlabel('Co-spike index')

if isempty(cosi)
    rho = nan; pval = nan;
else
    if do_rel
        [rho,pval] = corr(rel_increase,cosi,'Type','Spearman','rows','pairwise');
    else
        [rho,pval] = corr(abs_increase,cosi,'Type','Spearman','rows','pairwise');
    end
end
title(sprintf('%s rho = %1.2f p = %1.3f',name,rho,pval))
set(gca,'fontsize',20);

outfolder = [results_folder,'cospike/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end

print(gcf,[outfolder,name,'_cospike'],'-dpng');
close(gcf);


end