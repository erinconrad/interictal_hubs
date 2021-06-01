function [beta,pval] = rprep(thing,rate,change,surround,locs,spikey_idx,outfolder,...
    name,do_rel,ttext,labels,un_cos)

isnan_idx = isnan(thing);

thing = thing(spikey_idx & ~isnan_idx);
rate = rate(spikey_idx& ~isnan_idx,:);
locs = locs(spikey_idx& ~isnan_idx,:);
labels = labels(spikey_idx& ~isnan_idx);
un_cos = un_cos(spikey_idx& ~isnan_idx,spikey_idx& ~isnan_idx,:);

if isempty(surround)
    pre = nanmean(rate(:,1:change-1),2);
    post = nanmean(rate(:,change+1:end),2);
    un_cos_blocks = 1:size(un_cos,3);
else
    pre = nanmean(rate(:,change-surround:change-1),2);
    post = nanmean(rate(:,change+1:change+surround),2);
    un_cos_blocks = change-surround:change+surround;
    
end
if do_rel
    achange = (post-pre)./pre;
else
    achange = post-pre;
end

%uncosi = un_cos(:,:,un_cos_blocks)./sqrt(rate(:,un_cos_blocks))./(sqrt(rate(:,un_cos_blocks)))';

uncosi = nan(size(un_cos));
for i = 1:length(un_cos_blocks)
    uncosi(:,:,i) = un_cos(:,:,i)./sqrt(rate(:,i))./(sqrt(rate(:,i)))';
    
end

uncosi = nanmean(uncosi,3);
A = uncosi;
B = triu(A.',1) + tril(A);
B(B==inf) = 0;
B(isnan(B)) = 0;
% Correction because it is not perfectly symmetric (but it should be)


dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));

%fprintf('\n\nfirst order:\n');
[beta,pval] = spatial_autocorrelation(thing,achange,dmin,locs,B,1);

% simple corr
[rho_simp,pval_simp] = corr(thing,achange,'Type','Spearman','rows','pairwise');

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
title(sprintf(['simple corr: rho = %1.2f, p = %1.3f\n'...
    'regression coeff: %1.2f, p = %1.3f\nrho coeff:'...
    '%1.2f, p = %1.3f'],...
    rho_simp,pval_simp,beta(1),pval(1),beta(2),pval(2)));
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