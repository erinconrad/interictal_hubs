function [beta,pval] = rprep(thing,rate,change,surround,locs,spikey_idx,outfolder,...
    name,do_rel,ttext,labels)

isnan_idx = isnan(thing);

thing = thing(spikey_idx & ~isnan_idx);
rate = rate(spikey_idx& ~isnan_idx,:);
locs = locs(spikey_idx& ~isnan_idx,:);
labels = labels(spikey_idx& ~isnan_idx);

pre = nanmean(rate(:,change-surround:change-1),2);
post = nanmean(rate(:,change+1:change+surround),2);
if do_rel
    achange = (post-pre)./pre;
else
    achange = post-pre;
end



dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));

%fprintf('\n\nfirst order:\n');
[beta,pval] = spatial_autocorrelation(thing,achange,dmin,locs,1);

%{
fprintf('\n\nalternative:\n');
[beta,pval] = spatial_autocorrelation(thing,achange,dmin,locs,0);
%}


%% Find the 5 largest rate decrease electrodes
[sort_thing,I] = sort(achange);
big_thing_idx = I(1:5);
rate_big_thing = rate(big_thing_idx,:);
big_thing_labels = labels(big_thing_idx);

%
figure
set(gcf,'position',[100 300 1200 400])
t = tiledlayout(1,2);

nexttile
plot(thing,achange,'o','color',[1 1 1])
hold on
text(thing,achange,labels,'horizontalalignment','center','fontsize',20)
text(thing(big_thing_idx),achange(big_thing_idx),...
    labels(big_thing_idx),'horizontalalignment','center',...
    'fontsize',20,'color','r')
title(sprintf('regression coeff: %1.2f, p = %1.3f\nrho coeff: %1.2f, p = %1.3f',...
    beta(1),pval(1),beta(2),pval(2)));
ylabel('Spike rate change')
xlabel(ttext)
set(gca,'fontsize',20)

nexttile
plot(nanmean(rate_big_thing,1),'k')
hold on
plot([change change],ylim,'r--','linewidth',2)
xlabel('Block')
ylabel('Spike rate')
set(gca,'fontsize',20)

title(t,name,'fontsize',20);

print(gcf,[outfolder,name,'_',ttext,sprintf('_%d',surround)],'-dpng')
close(gcf)

%}


%{
T = table(locs(:,1),locs(:,2),locs(:,3),achange,rchange,thing,A,...
    'VariableNames',{'x','y','z','abs','rel','thing','A'});

writetable(T,[out_folder,'r_',name,'.csv'])

lm = fitlm(T,'abs~thing+A');
%}

end