function new = break_files_132(spikes)

new = spikes;
new_f = 0;
new_h = 0;


for f = 1:length(spikes.file)
    
    % advance new f with old file changes
    new_f = new_f + 1;
    subtract_time = 0;
 
    switch f
        case 1
            switch_time = nan;
        case 2
            switch_time = 86311.50;
        case 3
            switch_time = 336599.31;
            
    end
    
    if ~isnan(switch_time), did_switch = 0; else, did_switch = 1; end
    
    for h = 1:length(spikes.file(f).block)
        
        % advance new_h
        new_h = new_h + 1;
        
        start_time = spikes.file(f).block(h).run_times;
        
        % Figure out if I need to switch file and h
        if did_switch == 0 
            if start_time > switch_time
                new_f = new_f + 1; % advance new_f
                new_h = 1; % reset new_h
                did_switch = 1;
                subtract_time = start_time;
            end
        end
        
        new.file(new_f).block(new_h) = spikes.file(f).block(h);
        new.file(new_f).block(new_h).run_times = spikes.file(f).block(h).run_times - subtract_time;
        new_ch_labels = new.file(new_f).block(new_h).car_labels;
        
        % Make sure that, within a file, all channel labels are the same!
        if h == 1
            old_ch_labels = new_ch_labels; % set old channel labels to be those in the first block for the new file
        else
            if ~isequal(new_ch_labels,old_ch_labels)
                error('why');
            end
        end
        
    end
end

end