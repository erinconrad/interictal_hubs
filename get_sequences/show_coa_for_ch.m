function show_coa_for_ch(coa,ich,labels)

%% Accept strings
if ischar(ich)
    ich = find(strcmp(labels,ich));
elseif iscell(ich)
    ich = find(ismember(labels,ich));
end


%% If ndim coa > 2, mean across third dimension
if ndims(coa) == 3
    coa = nanmean(coa,3);
end

%% Get coa for ich
coa_ich = coa(ich,:)';

%% Sort
[sort_coa_ich,I] = sort(coa_ich,'descend');
sort_labels = labels(I);

%% Table
table(sort_labels,sort_coa_ich)

end