function [change,no_change_ever] = find_electrode_change_files(pt,p)


nfiles = length(pt(p).ieeg.file);
if nfiles == 1
    return
end

nchange = 0;

for f = 1:nfiles-1
    old_labels = clean_labels_2(pt(p).ieeg.file(f).chLabels(:,1));
    new_labels = clean_labels_2(pt(p).ieeg.file(f+1).chLabels(:,1));
    
    % A change file
    if ~isequal(old_labels,new_labels)
        
        lost_elecs = setdiff(old_labels,new_labels);
        added_elecs = setdiff(new_labels,old_labels);
        unchanged_elecs = intersect(new_labels,old_labels);
        
        if length(lost_elecs) + length(unchanged_elecs) ~= length(old_labels) || ...
                length(added_elecs) + length(unchanged_elecs) ~= length(new_labels)
            error('look');
        end
        
        nchange = nchange + 1;
        change(nchange).files = [f f+1];
        change(nchange).lost = lost_elecs;
        change(nchange).added = added_elecs;
        change(nchange).unchanged = unchanged_elecs;
        
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


end