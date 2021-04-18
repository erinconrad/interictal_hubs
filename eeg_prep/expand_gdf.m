function new_gdf = expand_gdf(gdf,chs_in_bipolar)

new_gdf = nan(size(gdf,2),3);
for i = 1:size(gdf,1)
    time = gdf(i,2);
    ch = gdf(i,1);
    new_chs = chs_in_bipolar(ch,:);
    
    index = (i-1)*2 + 1;
    new_gdf(index,:) = [new_chs(1) time i];
    new_gdf(index+1,:) = [new_chs(2) time i];
end

end