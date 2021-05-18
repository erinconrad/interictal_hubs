function new_distance(rel_increase,dist,spikey_idx,labels,name,results_folder)

if sum(spikey_idx) == 0
    fprintf('\nMissing spikey elecs, skipping\n');
    return
end

%% Restrict to spikey
rel_increase(~spikey_idx) = [];
dist(~spikey_idx) = [];
labels(~spikey_idx) = [];

thresh = [-0.5 1];
big = rel_increase > thresh(2);
small = rel_increase<thresh(1);
%[~,big] = max(rel_increase);
neither = ~big & ~small;



figure
set(gcf,'position',[440 512 1001 286])

%% Scatter
subplot(1,2,1)
text(rel_increase,dist,labels,...
    'horizontalalignment','center','fontsize',20);
hold on
text(rel_increase(big),dist(big),labels(big),...
    'horizontalalignment','center','fontsize',20,'color','g');
text(rel_increase(small),dist(small),labels(small),...
    'horizontalalignment','center','fontsize',20,'color','r');
xlim([min(rel_increase) max(rel_increase)])
ylim([0 max(dist)])
xlabel('Relative spike rate change')
ylabel('Distance from added')
set(gca,'fontsize',20)

%% Look at high increase (double or more) and big decrease (half or less)

subplot(1,2,2)
text(1+0.05*rand(sum(small),1),dist(small),labels(small),...
    'horizontalalignment','center','fontsize',20)
hold on
text(2+0.05*rand(sum(neither),1),dist(neither),labels(neither),...
    'horizontalalignment','center','fontsize',20)
text(3+0.05*rand(sum(big),1),dist(big),labels(big),...
    'horizontalalignment','center','fontsize',20)
xlim([0.5 3.5])
ylim([0 max(dist)])
xticks([1 2 3])
xticklabels({'Large decrease','Neither','Large increase'})
%xtickangle(45)
ylabel('Distance from added')

% ttest
group = ones(length(dist),1);
group(neither) = 2;
group(big) = 3;

pval = anova1(dist,group,'off');
title(sprintf('%s p = %1.3f',name,pval))
set(gca,'fontsize',20)

out_folder = [results_folder,'distance/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

print(gcf,[out_folder,name,'_distance'],'-dpng')
close(gcf)


end