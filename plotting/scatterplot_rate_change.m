function scatterplot_rate_change(rate_change,spikey_labels,spikes,p,name)

[sorted_change,I] = sort(rate_change);
sorted_labs = spikey_labels(I);

figure
plot(sorted_change)
xticks(1:length(sorted_change))
xticklabels(sorted_labs)
xtickangle(45)


while 1
    try
        [x,~] = ginput;
    catch
        break
    end
    chLab = sorted_labs{round(x(end))};
    fprintf('\nShowing spikes for %s ch %s\n',...
            name,chLab);
    plot_spikes_by_ch(p,chLab,[],[],spikes)
end
    

end