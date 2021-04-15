function sr = rl_stability(rl,nseq)

nchs = size(rl,1);
nhours = size(rl,2);

%% Median across segments
median_rl = nanmedian(rl,2);

sr = nan(nhours,1);

%% Get the spearman rank correlation between segment rl and overall rl
for h = 1:nhours

    r = corr(median_rl,rl(:,h),'Type','Spearman','Rows','pairwise');
    sr(h) = r;
end

if 0
figure
subplot(2,1,1)
plot(sr,'o')

subplot(2,1,2)
plot(nseq,'o')
end

end