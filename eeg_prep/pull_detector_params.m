function params = pull_detector_params(name,T)


%% Get the relevant row
row = find(strcmp(T.Patient,name));



if isempty(row)
    tmul = 17;
    absthresh = 50;
else
    
    tmul = T.Tmul(row);
    absthresh = T.Absthresh(row);

    %% Default values if not specified
    if isnan(tmul)
        tmul = 17;
    end
    if isnan(absthresh)
        absthresh = 50;
    end
end

params.tmul = tmul;
params.absthresh = absthresh;

end