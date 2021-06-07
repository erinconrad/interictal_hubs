function [gdf,n_removed]= remove_spikes_in_sz(gdf,sz_times)

% parameters
buffer = 0; % 0 seconds

% get times in one nice array
if isempty(sz_times)
    n_removed = 0;
    return
end

keep = ones(size(gdf,1),1);

% loop through spikes
for i = 1:size(gdf,1)
    spike_time = gdf(i,2);
    
    % if the spike is within 0 seconds of any start or end sz time
    if any(abs(spike_time-sz_times(:)) < buffer) || ...
            any(spike_time > sz_times(:,1) & spike_time < sz_times(:,2))
                    % if the spike is between the start and end (for long sz)

        keep(i) = 0; % throw it out
    end
end

keep = logical(keep);

gdf(~keep,:) = [];
n_removed = sum(~keep);


if sum(~keep) > 0
    fprintf('\nWarning, found %d spikes in seizures\n',sum(~keep))
end

end


