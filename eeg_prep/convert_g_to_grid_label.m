function out_labels = convert_g_to_grid_label(in_labels,name)

out_labels = cellfun(@(x) strrep(x,'RG','GRID'),in_labels,'UniformOutput',false);
out_labels = cellfun(@(x) strrep(x,'LG','GRID'),out_labels,'UniformOutput',false);

% Also convert AST to AS? What????
out_labels = cellfun(@(x) strrep(x,'AST','AS'),out_labels,'UniformOutput',false);
out_labels = cellfun(@(x) strrep(x,'MST','MS'),out_labels,'UniformOutput',false);

% Also convert LOF5 to LF5, ugh
out_labels = cellfun(@(x) strrep(x,'LOF5','LF5'),out_labels,'UniformOutput',false);

if strcmp(name,'HUP166')
   out_labels = cellfun(@(x) strrep(x,'LB10','LB1'),out_labels,'UniformOutput',false);
   out_labels = cellfun(@(x) strrep(x,'LB1','LB'),out_labels,'UniformOutput',false);
 
end

end