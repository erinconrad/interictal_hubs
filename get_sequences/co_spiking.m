function cos = co_spiking(gdf,fs,added_labels,chLabels)

%{
This measures, for each channel, how many times it has a spike at around
the same time one of the added channels has a spike
%}

t2 = 50*1e-3; % max time from preceding spike (15 ms in paper)

nchs = length(chLabels);
cos = zeros(nchs,1);

added_idx = find(ismember(chLabels,added_labels));

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
    
    % get spikes that are temporally close
    close = find(abs(time-gdf(:,2)) < t2);
    close(close == s) = []; % remove itself
    
    % get the spike chs
    close_spikes = gdf(close,1);
    
    % see if any are added channels
    if any(ismember(close_spikes,added_idx))
        
        % if so, add cos index
        cos(ch) = cos(ch) + 1;
    end
    
end

end