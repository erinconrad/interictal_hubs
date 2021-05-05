function raster_rate_chs(all_rate,last_block,block_dur,run_dur,change_block,...
    dist_info,unchanged,name,results_folder)

figure
    set(gcf,'position',[1 1 1400 800])
    times = 1:size(all_rate,2);
    times = times*block_dur;
    h = turn_nans_white(all_rate);
    set(h,'XData',[0:size(all_rate,2)*block_dur])
    xlim([0 size(all_rate,2)*block_dur])
    hold on
    %{
    for b = 1:length(last_block)
        plot([last_block(b) last_block(b)],ylim,'k--','linewidth',3)
    end
    %}
    cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',3);
    yticks(1:length(unchanged))
    title(sprintf('%s',name))
    ax = gca;
    set(gca,'yticklabel',unchanged)
    set(gca,'fontsize',10);


    ylabel('Hour')
    if 0
        while 1
            [x,y] = ginput;
            chLab = unchanged{round(y(end))};
            fidx = findices(round(x(end)));
            bidx = bindices(round(x(end)));
            fprintf('\nShowing spikes for %s ch %s file %d block %d\n',...
                pt_name,chLab,fidx,bidx);

            plot_spikes_by_ch(p,chLab,fidx,bidx)

        end
    end


    output_folder = [results_folder,'raster_rate/'];
    if exist(output_folder,'dir') == 0
        mkdir(output_folder)
    end

    print(gcf,[output_folder,name,'_raster'],'-dpng')

end