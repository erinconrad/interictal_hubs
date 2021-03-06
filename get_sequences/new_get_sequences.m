function [seq,rl,coa,num_seq,cos] = new_get_sequences(gdf,nchs,fs,...
    is_post,chLabels,added_labels,unchanged_labels)
    
t2 = 15*1e-3; % max time from preceding spike (15 ms in paper)
minSeqLength = 5;  % 5 in paper
%t2 = t2*fs;

ns = size(gdf,1);
coa = zeros(nchs,nchs);
num_seq = zeros(nchs,1);
cos = zeros(nchs,1);

%% Sort by time
times = gdf(:,2);
if ~isequal(times,sort(times))
    error('ruh-roh');
end

%% initialize seq cell
seq = {};

%% Initialize stuff for my while loop
curr_seq = [gdf(1,:),1];
next_index = 2;

while next_index <= ns
    
    next_time = gdf(next_index,2);
    last_time = curr_seq(end,2);
    
    % go to the next if the channel is already in the sequence
    if ismember(gdf(next_index,1),curr_seq(:,1))
        next_index = next_index+1;
        continue;
    end
    
    % if they're close enough in time
    if next_time - last_time < t2
        % append it to the sequence
        curr_seq = [curr_seq;[gdf(next_index,:) next_index]];
        
        % increase next_index
        next_index = next_index + 1;
        
    % if the next one is too far
    else
        
        
        % check if the current sequence is long enough
        if size(curr_seq,1) >= minSeqLength
            % Append it to the seq cell array
            seq = [seq;curr_seq];
            
            % Reset curr seq to be this new spike
            curr_seq = [gdf(next_index,:) next_index];
            
            % Increase next index
            next_index = next_index+1;
            
        % If it's too short, don't add it but increase index
        else
            % Reset curr seq to be this new spike
            curr_seq = [gdf(next_index,:) next_index];
            
            % Increase next index
            next_index = next_index+1;
            
        end
        
    end
    
end

nseq = length(seq);


rl = nan(nchs,nseq);

% Loop over sequences
for s = 1:nseq
    curr = seq{s};
    
    
    
    for ich = 1:nchs
        
        
        %% Get recruitment latency of every channel in each sequence
        if ismember(ich,curr(:,1))
            
            ich_index = find(curr(:,1) == ich);
            
            % difference between the time this electrode peaks and the
            % first peak in the sequence
            rl(ich,s) = curr(ich_index,2) - curr(1,2);
            
        end
        
        %{
        %% Get spike co-activation matrix
        % the co-activation for this segment between electrode i and j is the
        % number of sequences in which they co-occur
        for jch = 1:nchs
            % if i and j are in the sequence
            if ismember(ich,curr(:,1)) && ismember(jch,curr(:,1))
                coa(ich,jch) = coa(ich,jch) + 1; % increase coa for i,j and j,i by 1
                coa(jch,ich) = coa(jch,ich) + 1;
            end
        end
        %}

        
    end
    
    %
    %% Get spike coactivation matrix
    chs = curr(:,1);
    for i = 1:length(chs)
        num_seq(chs(i)) = num_seq(chs(i)) + 1;
        for j = 1:i-1
            ich = chs(i);
            jch = chs(j);
            coa(ich,jch) = coa(ich,jch) + 1;
            coa(jch,ich) = coa(jch,ich) + 1;
        end
    end
    %}
    
    if is_post
        
        % get channels idx corresponding to added labels
        added = find(ismember(chLabels,added_labels));
        
        % count if any added channel is co-spiking with the current channel
        chs = curr(:,1);
        for i = 1:length(chs)
            ich = chs(i);
            cos(ich) = cos(ich) + any(ismember(chs,added));
        end
        
        
    end
end

rl = nanmean(rl,2);


end

