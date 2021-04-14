function [which_chs,skip] = designate_chs(chLabels,clean_labs,clean_loc_labs,loc)

which_chs = (1:length(clean_labs))';

%% Match up channel labels
loc_indices = reconcile_labels(clean_labs,clean_loc_labs);

%% Skip electrodes marked "ignore" in localization struct
skip_ignore = loc_indices(loc.ignore == 1);

%% Skip electrodes marked outside of brain in localization struct
skip_extraaxial = loc_indices(loc.gm_wm == -1);

%% Skip ekg, rate, etc...
exclusion_text = {'ekg','EKG','ecg','ECG','rate','rr','RR'};

to_exclude = cellfun(@(x) contains(x,exclusion_text),chLabels,'UniformOutput',false);
to_exclude = cell2mat(to_exclude);

skip_ekg = which_chs(to_exclude);

%% Skip
skip.ignore = skip_ignore;
skip.extra_axial = skip_extraaxial;
skip.ekg = skip_ekg;
skip.all = unique([skip_ignore;skip_extraaxial;skip_ekg]);
which_chs(ismember(which_chs,[skip_ignore;skip_extraaxial;skip_ekg])) = [];


end

function loc_indices = reconcile_labels(clean_labs,clean_loc_labs)
% returns the index from the main channel labels corresponding to each of
% the localization labels

loc_indices = nan(length(clean_loc_labs),1);

for i = 1:length(clean_loc_labs)
    a = find(strcmp(clean_labs,clean_loc_labs{i}));
    if ~isempty(a)
        loc_indices(i) = a;
    end      
end

end