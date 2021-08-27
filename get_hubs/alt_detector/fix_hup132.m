function new_labels = fix_hup132(f,run_times,orig_labels,data_folder)

fname = [data_folder,'HUP132_fix.xlsx'];
start_time = run_times(1);

switch f
    case 1
        T = readtable(fname,'sheet','Original hookup');
    case 2
        if start_time < 86311.50
            T = readtable(fname,'sheet','Original hookup');
        else
            T = readtable(fname,'sheet','Replug 1');
        end        
    case 3
        if start_time < 336599.31
            T = readtable(fname,'sheet','Replug 2');
        else
            T = readtable(fname,'sheet','Revision');
        end
    
end
        
col_names = T.Properties.VariableNames;
orig_col = col_names{3};
new_col = col_names{4};

ieeg_labels = T.(orig_col);
corr_labels = T.(new_col);

%% Map ieeg_labels onto orig labels (should not have missing stuff)
% this is the mapping to say for each ieeg label, what is the corresponding
% original label index
[~,loc_i_to_o] = ismember(orig_labels,ieeg_labels);
if ~isequal(ieeg_labels(loc_i_to_o),orig_labels), error('why'); end

%% Get the corr_labels for these indices
new_labels = corr_labels(loc_i_to_o);

if 0
    table(orig_labels,ieeg_labels,corr_labels,new_labels)
end

end