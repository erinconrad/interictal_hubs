function [new,new_pc,new_ad] = break_files_132(spikes,pc,ad)

%new = spikes;
new_f = 0;
new_h = 0;
new.name = spikes.name;
new_pc.name = spikes.name;
new_ad.name = ad.name;

for f = 1:length(spikes.file)
    
    % advance new f and reset new h with old file changes
    new_f = new_f + 1;
    new_h = 0;
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
        
        if isempty(spikes.file(f).block(h).run_times), continue; end
        
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
        
        % fill in information from appropriate block into new struct
        new.file(new_f).block(new_h) = spikes.file(f).block(h);
        new.file(new_f).block(new_h).run_times = spikes.file(f).block(h).run_times - subtract_time;
        
        new_pc.file(new_f).block(new_h) = pc.file(f).block(h);
        new_pc.file(new_f).block(new_h).run_times = pc.file(f).block(h).run_times - subtract_time;
        
        new_ad.file(new_f).block(new_h) = ad.file(f).block(h);
        new_ad.file(new_f).block(new_h).run_times = new_ad.file(f).block(h).run_times - subtract_time;
        
        
        new_ch_labels = new.file(new_f).block(new_h).car_labels;
        
        % Make sure that, within a file, all channel labels are the same!
        if new_h == 1
            old_ch_labels = new_ch_labels; % set old channel labels to be those in the first block for the new file
        else
            if ~isequal(new_ch_labels,old_ch_labels)
                error('why');
            end
        end
        
    end
end

new = add_locs_hup132(new);

end