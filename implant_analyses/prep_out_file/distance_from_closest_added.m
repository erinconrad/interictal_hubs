function [dist,closest] = distance_from_closest_added(unchanged,added)

dist = nan(size(unchanged,1),1);
closest = nan(size(unchanged,1),1);

% Loop over unchanged electrodes
for i = 1:size(unchanged,1)
    
    added_dist = nan(size(added,1),1);
    
    % Loop over added electrodes
    for j = 1:size(added,1)
        
        % Calculate distance between the two
        added_dist(j) = vecnorm(unchanged(i,:) - added(j,:));
        
    end
    
    % Find the minimum dist
    [dist(i),closest(i)] = min(added_dist);
    
end

end