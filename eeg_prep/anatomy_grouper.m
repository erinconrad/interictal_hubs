function [lat,loc] = anatomy_grouper(anatomy)

loc = cell(length(anatomy),1);
lat = cell(length(anatomy),1);

for ich = 1:length(anatomy)
    curr = anatomy{ich};
    
    if isempty(curr)
        lat{ich} = 'Unspecified';
        loc{ich} = 'Unspecified';
        continue
    end
    
    
    %% Get laterality
    if contains(curr,'right','IgnoreCase',true) || strcmp(curr(1),'r') || strcmp(curr(1),'R')
        lat{ich} = 'Right';
    end
    
    if contains(curr,'left','IgnoreCase',true) || strcmp(curr(1),'l') || strcmp(curr(1),'L')
        if strcmp(lat{ich},'R')
            error('check laterality');
        else
            lat{ich} = 'Left';
        end
    end
    
    %% Get localization
    if contains(curr,'white','IgnoreCase',true)
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'white matter';
    end
    
    if contains(curr,'amy','IgnoreCase',true) || ...
            contains(curr,'hipp','IgnoreCase',true) || ...
            contains(curr,'temp','IgnoreCase',true) && ...
            (contains(curr,'med','IgnoreCase',true) || ...
            contains(curr,'mes','IgnoreCase',true))
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'mesial temporal';
    end

    if contains(curr,'temp','IgnoreCase',true) && ...
            (~contains(curr,'med','IgnoreCase',true) && ...
            ~contains(curr,'mes','IgnoreCase',true))
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'temporal neocortical';
    end
    
    %{
    if contains(curr,'front','IgnoreCase',true)
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'frontal';
    end
    %}
    
    % if don't find it, fill it up
    if isempty(loc{ich})
        loc{ich} = 'other';
    end
    
    
end

end