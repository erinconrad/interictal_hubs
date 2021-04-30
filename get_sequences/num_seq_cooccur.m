function perc = num_seq_cooccur(coa,nseq,ich,jch,chLabels)

%% Accept strings
if ischar(ich)
    ich = find(strcmp(chLabels,ich));
elseif iscell(ich)
    ich = find(ismember(chLabels,ich));
end

if ischar(jch)
    jch = find(strcmp(chLabels,jch));
elseif iscell(jch)
    jch = find(ismember(chLabels,jch));
end

%% If ndim coa > 2, mean across third dimension
if ndims(coa) == 3
    coa = nanmean(coa,3);
end

if ndims(nseq) == 2
    nseq = nanmean(nseq,2);
end

perc = zeros(length(ich),1);

for i = 1:length(ich)
    
    curr_ich = ich(i);
    
    % find the highest co-occuring channel amongst jchs, and get how often this
    % happens
    num_coa = max(coa(curr_ich,jch));

    % number of total sequences with ich
    num_seq = nseq(curr_ich);

    % what percent of sequences involving ich also involve its maximally
    % co-occurring ch within this list
    perc(i) = num_coa/num_seq;
end
    
    
        
end

