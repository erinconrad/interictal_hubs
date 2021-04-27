function T = find_szs_in_annotations(pt)

poss_sz_text = {};
poss_sz_start = {};
files = [];
names = [];
for p = 1:length(pt)
    
    if isempty(pt(p).ieeg), continue; end
    
    for f = 1:length(pt(p).ieeg.file)
        if ~isfield(pt(p).ieeg.file(f),'ann'), continue; end
        if strcmp(pt(p).ieeg.file(f).ann,'empty'), continue; end
        n_anns = length(pt(p).ieeg.file(f).ann);
        for a = 1:n_anns
            
            n_events = length(pt(p).ieeg.file(f).ann(a).event);
            for i = 1:n_events


                description = pt(p).ieeg.file(f).ann(a).event(i).description;

                % search for seizure-y strings
                if contains(description,'seizure','IgnoreCase',true) || ...
                        contains(description,'sz','IgnoreCase',true) || ...
                        contains(description,'onset','IgnoreCase',true) || ...
                        contains(description,'UEO','IgnoreCase',true) || ...
                        contains(description,'EEC','IgnoreCase',true) 

                    poss_sz_text = [poss_sz_text;description];
                    poss_sz_start = [poss_sz_start;...
                        sprintf('%d',floor(pt(p).ieeg.file(f).ann(a).event(i).start))];
                    files = [files;f];
                    names = [names;pt(p).name];

                end
            end

        end



    end
end


T = table(names,files,...
    poss_sz_start,...
    poss_sz_text);

end