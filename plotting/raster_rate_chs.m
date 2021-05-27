function raster_rate_chs(all_rate,block_dur,change_block,...
    unchanged,name,results_folder,findices,bindices,p,spikes,...
    spikey_idx,surround)

figure
    set(gcf,'position',[1 1 1400 800])
    times = 1:size(all_rate,2);
    times = times*block_dur;
    h = turn_nans_white(all_rate);
    set(h,'XData',[0:size(all_rate,2)*block_dur])
    xlim([0 size(all_rate,2)*block_dur])
    hold on
    cp = plot([change_block*block_dur change_block*block_dur],ylim,'r--','linewidth',3);
    yticks(1:length(unchanged))
    title(sprintf('%s',name))
    ax = gca;
    set(gca,'yticklabel',unchanged)
    set(gca,'fontsize',10);


    xlabel('Hour')
    
    
    %% Permanova analysis
    if 1
        relative = 1;
        if relative == 0
            rel_text = 'Pre-post difference';
        else
            rel_text = 'Pre-post relative difference';
        end
        pval = ec_permanova(all_rate,spikey_idx,change_block,surround,relative);
        title(sprintf('%s %s p = %1.3f',name,rel_text,pval))

        output_folder = [results_folder,'raster_rate/'];
        if exist(output_folder,'dir') == 0
            mkdir(output_folder)
        end

    
        if relative
            print(gcf,[output_folder,name,'_rel'],'-dpng')
        else
            print(gcf,[output_folder,name,'_abs'],'-dpng')
        end
        close(gcf)
    end
    
    
    %% Pick a desired ch and block and show spikes
    if 0
        while 1
            try
                [x,y] = ginput;
            catch
                return
            end
            chLab = unchanged{round(y(end))};
            fidx = findices(round(x(end)/block_dur));
            bidx = bindices(round(x(end)/block_dur));
            fprintf('\nShowing spikes for %s ch %s file %d block %d\n',...
                name,chLab,fidx,bidx);

            plot_spikes_by_ch(p,chLab,fidx,bidx,spikes)

        end
    end


    

end