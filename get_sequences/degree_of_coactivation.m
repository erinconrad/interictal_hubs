function [ich_deg,jch_deg] = degree_of_coactivation(coa,ich,jch)

ich_deg = sum(coa(ich,jch))/sum(coa(ich,:));
jch_deg = sum(coa(ich,jch))/sum(coa(jch,:));

end