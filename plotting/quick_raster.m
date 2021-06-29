function quick_raster(out,p)

if ischar(p)
    for i = 1:length(out)
        if strcmp(p,out(i).name)
            p = i;
            break
        end
    end
end
        

rate = out(p).rate;
chLabels = out(p).unchanged_labels;
cblock = out(p).change_block;

ekg = identify_ekg_scalp(chLabels);

rate = rate(~ekg,:);
chLabels = chLabels(~ekg);


figure
turn_nans_white(rate)
hold on
plot([cblock cblock],ylim,'r--','linewidth',4)
yticks(1:length(chLabels))
yticklabels(chLabels)

%{
while 1
    try
        [x,y] = ginput;
    catch
        break
    end
    chLab = chLabels{round(y(end))};
    fidx = f;
    bidx = bindices(round(x(end)/dur));
    fprintf('\nShowing spikes for %s ch %s file %d block %d (hour %d)\n',...
        name,chLab,fidx,bidx,round(bidx*dur));

    plot_spikes_by_ch(p,chLab,fidx,bidx,spikes)

end
%}


end