function do_nmf(coa)



nchs = size(coa,1);
nblocks = size(coa,3);

 %% Flatten coa matrix
A = nan(nchs*(nchs-1)/2,nblocks);
for j = 1:size(coa,3)
    A(:,j) = wrap_or_unwrap_adjacency(coa(:,:,j));
end

A(isnan(A)) = 0;

%% Do NMF
k = 2;
[W,H] = nnmf(A,k);

%% Unpack subgraphs
subgraphs = zeros(nchs,nchs,k);
for j = 1:k
    subgraphs(:,:,j) = wrap_or_unwrap_adjacency(W(:,j));
end

%% Plot subgraphs
if 1
figure
for j = 1:k
    subplot(1,k,j)
    imagesc(subgraphs(:,:,j))
    title(sprintf('Subgraph %d',j))
end
end

%% Plot change in subgraphs over time
if 1
figure
lp = zeros(k,1);
lpn = cell(k,1);
for j = 1:k
    lp(j) = plot(H(j,:));
    lpn{j} = sprintf('Subgraph %d',j);
    hold on
end
legend(lp,lpn)

end

end