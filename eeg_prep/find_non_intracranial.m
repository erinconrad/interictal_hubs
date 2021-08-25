function ekg = find_non_intracranial(labels)

ekg = zeros(length(labels),1);

for i = 1:length(labels)
    
    if contains(labels(i),'ekg','ignorecase',true)
        ekg(i) = 1;
    end
    
    if contains(labels(i),'ecg','ignorecase',true)
        ekg(i) = 1;
    end
    
    if contains(labels(i),'rate','ignorecase',true)
        ekg(i) = 1;
    end
    
    if contains(labels(i),'rr','ignorecase',true)
        ekg(i) = 1;
    end
    
    if strcmp(labels(i),'C3') || strcmp(labels(i),'C4') || ...
            strcmp(labels(i),'CZ') || ...
            strcmp(labels(i),'F8') || ...
            strcmp(labels(i),'F4') || ...
            strcmp(labels(i),'FP2') || ...
            strcmp(labels(i),'FZ') || ...
            strcmp(labels(i),'LOC') || ...
            strcmp(labels(i),'T4') || ...
            strcmp(labels(i),'C6') || ...
            strcmp(labels(i),'ROC') || ...
            strcmp(labels(i),'P4') || ...
            strcmp(labels(i),'T6')

        ekg(i) = 1;
    end
    
    % fix for things that could be either scalp or ieeg
    if strcmp(labels(i),'O2')
        if sum(strcmp(labels,'O1')) == 0 % if hemiscalp, should not have odd; if ieeg, should have O1
            ekg(i) = 1;
        end
    end

    
end

ekg = logical(ekg);

end