function times = all_sz_times_in_file(pt,p,f)

if ~isfield(pt(p).ieeg.file(f),'sz')
    times = [];
    return;
end

sz = pt(p).ieeg.file(f).sz;
if isempty(sz)
    times = [];
    return
end

times = nan(length(sz),2);

for s = 1:length(sz)
    eec = sz(s).eec;
    ueo = sz(s).ueo;
    send = sz(s).end;
    
    % preferentially use eec
    if isnan(eec)
        times(s,:) = [ueo send];
    else
        times(s,:) = [eec send];
    end
    
end

end