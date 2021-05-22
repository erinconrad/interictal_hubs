function rprep(thing,rate,change,surround,locs,spikey_idx,results_folder,...
    name,do_rel)

thing = thing(spikey_idx);
rate = rate(spikey_idx,:);
locs = locs(spikey_idx,:);

pre = nanmean(rate(:,change-surround:change-1),2);
post = nanmean(rate(:,change+1:change+surround),2);
if do_rel
    achange = (post-pre)./pre;
else
    achange = post-pre;
end

out_folder = [results_folder,'for_r/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder);
end



dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));

spatial_autocorrelation(thing,achange,dmin,locs);


%{
T = table(locs(:,1),locs(:,2),locs(:,3),achange,rchange,thing,A,...
    'VariableNames',{'x','y','z','abs','rel','thing','A'});

writetable(T,[out_folder,'r_',name,'.csv'])

lm = fitlm(T,'abs~thing+A');
%}

end