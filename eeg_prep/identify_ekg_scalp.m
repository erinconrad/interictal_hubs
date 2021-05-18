function ekg = identify_ekg_scalp(labels)

ekg = zeros(length(labels),1);

for i = 1:length(labels)
    
    if contains(labels(i),'ekg','ignorecase',true)
        ekg(i) = 1;
    end
    
    if contains(labels(i),'ecg','ignorecase',true)
        ekg(i) = 1;
    end
    
    if strcmp(labels(i),'C3') || strcmp(labels(i),'C4') || ...
            strcmp(labels(i),'CZ') || ...
            strcmp(labels(i),'F8') || ...
            strcmp(labels(i),'FZ') || ...
            strcmp(labels(i),'LOC') || ...
            strcmp(labels(i),'O2') || ...
            strcmp(labels(i),'T4') || ...
            strcmp(labels(i),'C6') || ...
            strcmp(labels(i),'ROC')
        ekg(i) = 1;
    end
    
end

ekg = logical(ekg);

end