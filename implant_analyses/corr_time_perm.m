function tstat = corr_time_perm(thing,rate,change,surround,spikey_idx,do_abs,...
    ttext,outfolder,labels,name,do_simple)

%{
The general idea is that under the null distribution, the change in any
statistic from pre-to-post revision is not any larger than the change in
any statistic from pre-to-post any randomly chosen time.

In this particular case, the correlation between the spike rate change and
the "thing", where thing may be distance from the nearest added electrode
or co-spike index or whatever, is not expected to be any higher (under the
null) when looking around the revision time compared to any other time.

A significant result means that implantation has a particular tendency to
affect the spike rate preferentially for special electrodes, and that this
effect is stronger than the effect from other random times (sleep state,
etc.)

A non-significant result allows for the possibility that implantation may
preferentially increase spike rates in special electrodes, but implies that
it does not do so more than other times do.
%}

%% Restrict to spikey
if 1
thing = thing(spikey_idx);
rate = rate(spikey_idx,:);
labels = labels(spikey_idx);
end

if sum(spikey_idx) == 0
    fprintf('\nSkipping %s due to low spike rates\n',name);
    return
end

%% Get true rate change
pre = nanmean(rate(:,change-surround:change-1),2);
post = nanmean(rate(:,change+1:change+surround),2);

%% Decide whether to do absolute change or relative change
if do_abs
    rchange = post-pre;
    atext = 'abs';
else
    rchange = (post-pre)./pre;
    atext = 'rel';
end

if do_simple
    stext = 'simple';
else
    stext = '';
end

if sum(~isnan(rchange)) == 0
    fprintf('\nSkipping %s as rchange is all nans\n',name);
    return
end

%% Correlate change with thing
% thing could be distance from nearest added electrode, or co-spike index
[true_rho,pval_simple] = corr(rchange,thing,'Type','Spearman','rows','pairwise');
n = sum(~isnan(rchange));
tstat = true_rho * sqrt(n-2)/sqrt(1-true_rho^2);

%% Start permutation stuff
nb = 1e3;
perm_rho = nan(nb,1);
nblocks = size(rate,2);

% Loop over permutations
for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([change,nblocks-surround]);
    
    % recalculate change around this time
    fpre = nanmean(rate(:,fchange-surround:fchange-1),2);
    fpost = nanmean(rate(:,fchange+1:fchange+surround),2);
    if do_abs
        fchange = fpost-fpre;
    else
        fchange = (fpost-fpre)./fpre;
    end
    
    % correlate thing with fake change
    fcorr = corr(fchange,thing,'Type','Spearman','rows','pairwise');
    perm_rho(ib) = fcorr;
    
end

%% Sort the permutation results
perm_rho = sort(perm_rho);

%% Determine number with as or more extreme a correlation
num_as_sig = sum(abs(perm_rho) >= abs(true_rho)); % those more positive or more negative

pval = (num_as_sig+1)/(nb+1);

%% plot
if 0
    figure
    plot(perm_rho,'o')
    hold on
    plot(xlim,[true_rho true_rho])
    title(sprintf('p = %1.3f',pval))
end

figure
set(gcf,'position',[440 512 1001 286])
text(rchange,thing,labels,...
    'horizontalalignment','center','fontsize',20);
hold on
xlim([min(rchange) max(rchange)])
ylim([min(thing) max(thing)])
xlabel('Spike rate change')
ylabel(ttext)
if do_simple
    title(sprintf('%s rho = %1.2f (p = %1.3f)',name,true_rho,pval_simple))
else
    title(sprintf('%s rho = %1.2f (p = %1.3f by permutation test)',name,true_rho,pval))
end
set(gca,'fontsize',20)

print(gcf,[outfolder,name,'_',ttext,'_',sprintf('%d',surround)],'-dpng')
close(gcf)



end