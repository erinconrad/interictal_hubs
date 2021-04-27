function all_missed = check_for_missing_szs(pt,p)

all_missed = 0;

if isempty(pt(p).seizure_info)
    all_missed = 1;
    return;
end


num_missing = 0;
for i = 1:length(pt(p).seizure_info.sz)
    if strcmp(pt(p).seizure_info.sz(i).UEO,'missing')
        num_missing = num_missing + 1;
    end
    
end

if num_missing == length(pt(p).seizure_info.sz)
    all_missed = 1;
    return;
end

end