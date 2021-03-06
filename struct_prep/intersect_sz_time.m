function [int,all_inside] = intersect_sz_time(run_times,szs)

int = 0;
all_inside = 0;

for s = 1:length(szs)
    eec = szs(s).EEC;
    end_time = szs(s).End;
    
    if strcmp(eec,'missing')
        eec = szs(s).UEO;
    end
    
    if strcmp(eec,'missing')
        continue;
    end
    
    % fix for silliness
    if isempty(end_time) || end_time < eec
        end_time = eec+100;
    end
    
    % is the run contained in the seizure
    if run_times(1) >= eec && run_times(2) <= end_time
        int = 1;
        all_inside = 1;
    % does the run start during the seizure
    elseif run_times(1) >= eec && run_times(1) <= end_time
        int = 1;
    % does the run end during the seizure
    elseif run_times(2) >= eec && run_times(2) <= end_time
        int = 1;
    % does the run contain the seizure
    elseif run_times(1) <= eec && run_times(2) >= end_time
        int = 1;    
    end
end

end