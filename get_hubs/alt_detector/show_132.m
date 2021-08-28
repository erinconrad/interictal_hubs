function show_132(out)

rate = out.rate;
labels = out.unchanged_labels;
findices = out.findices;
fchange = find(diff(findices)~=0);

turn_nans_white(rate);
hold on
yticks(1:length(labels))
yticklabels(labels);
for i = 1:length(fchange)
    plot([fchange(i) fchange(i)],ylim,'r--','linewidth',3)
end

plot([out.change_block out.change_block],ylim,'y','linewidth',3)

end