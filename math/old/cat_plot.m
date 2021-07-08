function catlist = cat_plot(list1,list2)

if length(list1) ~= length(list2)
    error('lengths must be equal');
end

n = length(list1);
catlist = nan(n,1);


% turn nans into negative infs
list1(isnan(list1)) = -inf;
list2(isnan(list2)) = -inf;

% Sort them in descending order
[~,list1o] = sort(list1,'descend');
[~,list2o] = sort(list2,'descend');

for i = 1:n
    catlist(i) = length(intersect(list1o(1:i),list2o(1:i)))/i;
end


end