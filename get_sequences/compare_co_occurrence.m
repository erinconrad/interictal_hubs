function compare_co_occurrence(coa,unchanged,added,labels)

if iscell(unchanged)
    unchanged = find(ismember(labels,unchanged));
end

if iscell(added)
    added = find(ismember(labels,added));
end

%% If ndim coa > 2, mean across third dimension
if ndims(coa) == 3
    coa = nanmean(coa,3);
end



%% Identify electrodes with a large increase in spike rate


%% For each unchanged electrode, identify its maximally co-occurring added electrode
nunchanged = length(unchanged);
all_max_add = nan(nunchanged,1);
all_added_rank = nan(nunchanged,1);

for i = 1:nunchanged
    ich = unchanged(i);
    
    % co-activation array between ich and all added electrode
    coa_added = coa(ich,added);
    
    % this is the index of the added electrode with max co-occurrence for
    % this unchanged electrode
    [~,max_add] = max(coa_added);
    max_add = added(max_add);
    all_max_add(i) = max_add;
    
    % Make a list comprise of other unchanged electrodes and this maximally
    % co-occurring added electrode
    un_plus = [max_add;unchanged]; % 
    
    % Get the coas of this list
    coa_un_plus = coa(ich,un_plus);
    
    if sum(coa_un_plus) == 0
        added_rank = nan;
    else
    
        % sort in descending order
        [~,I] = sort(coa_un_plus,'descend');

        % Get the ranking of the added electrode (I expect this would be higher
        % for the sentinel electrode) - figure out how to deal with ties!!!!!
        added_rank = find(I == 1);
        
    end
    
    all_added_rank(i) = added_rank;
    
end

[sorted_ranks,I] = sort(all_added_rank);

T = table(labels(unchanged(I)),labels((all_max_add(I))),sorted_ranks);


end