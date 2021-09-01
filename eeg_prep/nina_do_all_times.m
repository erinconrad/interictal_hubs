function pt = nina_do_all_times(pt)

for p = 1:length(pt)
    if isempty(pt(p).ieeg),continue; end
    for f = 1:length(pt(p).ieeg.file)
        for b = 1:length(pt(p).ieeg.file(f).block)
            pt(p).ieeg.file(f).block(b).run = ...
                [pt(p).ieeg.file(f).block(b).start,....
                pt(p).ieeg.file(f).block(b).end];
        end
    end

end

end