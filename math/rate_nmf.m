function [W,H] = rate_nmf(rate,k)

rate(isnan(rate)) = 0;
A = rate;

fnorm = nan(k,1);

%%
for j = 1:k

    %% Do NMF
    [W,H] = nnmf(A,j);
    
    %% Compute F norm
    fnorm(j) = norm(A-W*H,'fro');
    
end

end