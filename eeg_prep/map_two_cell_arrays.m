function map_two_cell_arrays(A,B)

% Find idx such that 

[tf,loc]=ismember(A,B);
idx=[1:length(A)];
idx=idx(tf);
idx=idx(loc(tf));

if ~isequal(A(idx),B), error('why'); end

end