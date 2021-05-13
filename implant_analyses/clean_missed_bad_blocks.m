function spikes = clean_missed_bad_blocks(spikes)

%% Parameters
num_to_search = [-3 3];
max_nan = 3;

% Loop across files
for f = 1:length(spikes.file)
    
    nblocks = length(spikes.file(f).block);
    for ib = -num_to_search(1) + 1 : nblocks - num_to_search(2)
        
        num_total_skip = 0;
        for s = ib + num_to_search(1) : ib + num_to_search(2)
            if spikes.file(f).block(s).run_skip == 1
                num_total_skip = num_total_skip + 1;
            end
        end
        
        % if many surrounding are bad
        if num_total_skip > max_nan

            % call this one bad too
            spikes.file(f).block(ib).run_skip = 1;
            
            % don't do channel thing
            continue

        end
        
        nchs = length(spikes.file(f).block(ib).chLabels);
        for ich = 1:nchs
            
            num_nan = 0;
            
            % Look at the surrounding few blocks
            for s = ib + num_to_search(1) : ib + num_to_search(2)
                
                % add if it's bad
                if ismember(ich,spikes.file(f).block(s).bad)
                    num_nan = num_nan + 1;
                end
                
            
            end
            
            % if many surrounding are bad
            if num_nan > max_nan
                
                % call this one bad too
                spikes.file(f).block(ib).bad = unique([spikes.file(f).block(ib).bad;ich]);
                
            end
            
        end
        
        
    end
    
end


end