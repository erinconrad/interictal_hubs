function [cosi,added_rate] = alt_co_spiking(gdf,added_labels,unchanged_labels,chLabels)

t2 = 50*1e-3; % max time from preceding spike (15 ms in paper)

% cos(i,j) will be the number of times in which there is a co-spiking event
% between added channel i and unchanged channel j
cos = zeros(length(added_labels),length(unchanged_labels));
added_rate = zeros(length(added_labels),1);

added_chs = find(ismember(chLabels,added_labels));
unchanged_chs= find(ismember(chLabels,unchanged_labels));

%% Sort by time
times = gdf(:,2);
if ~isequal(times,sort(times))
    error('ruh-roh');
end

% Loop over spikes
for s = 1:size(gdf,1)
    
    % get time
    time = gdf(s,2);
    
    % get ch
    ch = gdf(s,1);
    
    % skip if it's not an added channel
    if ~ismember(ch,added_chs)
        continue
    end
    
    added_idx = find(ismember(added_chs,ch));
    added_rate(added_idx) = added_rate(added_idx) + 1;
    
    % get spikes that are temporally close
    close = find(abs(time-gdf(:,2)) < t2);
    close(close == s) = []; % remove itself
    
    % get the spike chs
    close_spikes = gdf(close,1);
    close_spikes = unique(close_spikes);
    
    % loop over the close spikes
    for c = 1:length(close_spikes)
        
        cch = close_spikes(c);
        
        % if it's an unchanged channel
        if ismember(cch,unchanged_chs)
            
            unchanged_idx = find(ismember(unchanged_chs,cch));
            
            % add it
            cos(added_idx,unchanged_idx) = cos(added_idx,unchanged_idx) + 1;
        end
    end
    
end

cosi = cos./added_rate;

end