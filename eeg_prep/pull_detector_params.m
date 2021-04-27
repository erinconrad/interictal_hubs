function params = pull_detector_params(name,T)


%% Get the relevant row
row = find(strcmp(T.Patient,name));

tmul = T.Tmul(row);
absthresh = T.Absthresh(row);

%% Default values if not specified
if isnan(tmul) || isempty(row)
    tmul = 17;
end
if isnan(absthresh) || isempty(row)
    absthresh = 50;
end

params.tmul = tmul;
params.absthresh = absthresh;

end