function perc_diff = comp_sequences(seq,chs,comp_chs,chLabels)

%% Accept strings
if ischar(chs)
    chs = find(strcmp(chLabels,chs));
elseif iscell(chs)
    chs = find(ismember(chLabels,chs));
end

if ischar(comp_chs)
    comp_chs = find(strcmp(chLabels,comp_chs));
elseif iscell(comp_chs)
    comp_chs = find(ismember(chLabels,comp_chs));
end

nchs = length(chs);
ncomp = length(comp_chs);

nseq = zeros(nchs,1);
nshare = zeros(nchs,nchs);
nbefore = zeros(nchs,nchs);
nafter = zeros(nchs,nchs);

for s = 1:length(seq)
    curr = seq{s};
    for i = 1:size(curr,1)
        ch = curr(i,1);
        nseq(ch) = nseq(ch) + 1;
        
        % Look at preceding spikes in sequence
        for j = 1:i-1
            
            % if it's a comparison electrode
            jch = curr(j,1);
            if ismember(jch,comp_chs)
                
                nshare(ch,jch) = nshare(ch,jch) + 1;
                
                % IMPORTANT - make sure the time is truly earlier (rather
                % than equal)
                if curr(j,2) < curr(i,2)
                    nbefore(ch,jch) = nbefore(ch,jch) + 1;
                end
                
            end
            
        end
        
        % now Look at later spikes in sequence
        for j = i+1:size(curr,1)
            
            % if it's a comparison electrode
            jch = curr(j,1);
            if ismember(jch,comp_chs)
                nshare(ch,jch) = nshare(ch,jch) + 1;
                
                % IMPORTANT - make sure the time is truly later (rather
                % than equal)
                if curr(j,2) > curr(i,2)
                    nafter(ch,jch) = nafter(ch,jch) + 1;
                end
                
            end
            
        end
    end
    
end

perc_before = nbefore./nshare;
perc_after = nafter./nshare;

perc_diff = perc_before - perc_after;

end