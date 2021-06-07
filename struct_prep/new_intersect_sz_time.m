function [int,all_inside] = new_intersect_sz_time(run_times,szs)

int = 0;
all_inside = 0;

for s = 1:size(szs,1)
    eec = szs(s,1);
    end_time = szs(s,2);
    
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