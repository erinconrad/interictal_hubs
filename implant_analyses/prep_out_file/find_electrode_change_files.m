function [change,no_change_ever] = find_electrode_change_files(pt,p,only_depth)


nfiles = length(pt(p).ieeg.file);
if nfiles == 1
    return
end

nchange = 0;
name = pt(p).name;

for f = 1:nfiles-1
    old_labels = clean_labels_2(pt(p).ieeg.file(f).chLabels(:,1));
    new_labels = clean_labels_2(pt(p).ieeg.file(f+1).chLabels(:,1));
    
    % A change file
    if ~isequal(old_labels,new_labels)
        
        lost_elecs = setdiff(old_labels,new_labels);
        added_elecs_f = setdiff(new_labels,old_labels);
        unchanged_elecs = intersect(new_labels,old_labels);
        
        if length(lost_elecs) + length(unchanged_elecs) ~= length(old_labels) || ...
                length(added_elecs_f) + length(unchanged_elecs) ~= length(new_labels)
            error('look');
        end
        
        nchange = nchange + 1;
        change(nchange).files = [f f+1];
        change(nchange).old_lost = lost_elecs;
        change(nchange).old_added = added_elecs_f;
        change(nchange).old_unchanged = unchanged_elecs;
        
        [new_elecs,new_depths,really_unchanged] = added_elecs(name,new_labels,unchanged_elecs);
        new_subdural = new_elecs(~ismember(new_elecs,new_depths));
        if only_depth
            change(nchange).added = new_depths;
        else
            change(nchange).added = new_elecs;
        end
        change(nchange).added_depth = new_depths;
        change(nchange).added_subdural = new_subdural;
        change(nchange).added_all = new_depths;
        change(nchange).unchanged = really_unchanged;
        
        
    end
end

all_unchanged = cell(length(change),1);
for c = 1:length(change)
    all_unchanged{c} = change(c).unchanged;
    
end

no_change_ever = all_unchanged{1};
for c = 2:length(change)
    no_change_ever = intersect(no_change_ever,all_unchanged{c});
end

if 0
    for c = 1:length(change)
        fprintf('\nFor %s, unchanged elecs are:\n',name);
        change(c).unchanged
    
        fprintf('\nAdded elecs are:\n');
        change(c).added
    end
    
end


end