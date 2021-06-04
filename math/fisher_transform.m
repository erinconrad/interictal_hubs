function [bizarro_z,z_score,pval] = fisher_transform(r,n)

bizarro_z = atanh(r);
ste = sqrt(n-3);
z_score = bizarro_z*ste;
pval = 2*normcdf(-abs(z_score));

end