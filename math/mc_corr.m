function [true_rho,pval,mc_rho] = mc_corr(rate,ns,predictor,cblock,surround,nb,which_resp,corr_type,only_pre)

%{
This is the main statistical test for the correlation analysis of the
implant effect paper. Essentially, I pick a random pseudo-revision time, I
recalculate the correlation between the rate change and the distance, and I
see how many pseudo-revision times yield as or more extreme a correlation
as the true correlation.
%}

nblocks = size(rate,2);

%% Get true corr
% Identify pre and post times
[pre,post] = get_surround_times(rate,cblock,surround);

% Calculate change in rate
switch which_resp
    case 'rel_rate'
        resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
    case 'ns_rel'
        resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./nanmean(ns(:,pre),2);
end

true_rho = corr(resp,predictor,'Type',corr_type,'rows','pairwise');

mc_rho = nan(nb,1);
for ib = 1:nb
    
    % I wrap this in a while loop because some choices of times result in
    % nan correlations. I don't want to accept these because I can't
    % compare them with my true correlation. And so if I get one I choose a
    % new time.
    while 1
        
        % Make a fake change time
        if only_pre
            fchange = randi([surround+1,cblock-1]);
        else
            fchange = randi([surround+1,nblocks-surround]);
        end

        % recalculate change around this time
        [fpre,fpost] = get_surround_times(rate,fchange,surround);
    
        if length(fpre) == 1 || length(fpost) == 1 
            continue % bad time, try another
        end
    
    
        switch which_resp
            case 'rel_rate'
                resp = (nanmean(rate(:,fpost),2) - nanmean(rate(:,fpre),2))./nanmean(rate(:,fpre),2);
            case 'ns_rel'
                resp = (nanmean(ns(:,fpost),2) - nanmean(ns(:,fpre),2))./nanmean(ns(:,fpre),2);
        end

        mc_rho_temp = corr(resp,predictor,'Type',corr_type,'rows','pairwise');
        
        if isnan(mc_rho_temp)
            continue % bad time, try another
        else
            break % if not a nan, cool, accept it
        end
    
    end
    
    mc_rho(ib) = mc_rho_temp;
end


%% Determine number with as or more extreme a correlation
% two sided
num_as_sig = sum(abs(mc_rho) >= abs(true_rho)); % those more positive or more negative

pval = (num_as_sig+1)/(nb+1);

% I should not have any nans
if sum(isnan(mc_rho))>0
    error('why');
end

if isnan(true_rho)
    error('why');
end

[sortedmc,I] = sort(mc_rho);


if 0
    figure
    plot(sortedmc)
    hold on
    if true_rho < 0
        plot(find(sortedmc<=true_rho),sortedmc(sortedmc<=true_rho),'r')
    else
        plot(find(sortedmc>=true_rho),sortedmc(sortedmc>=true_rho),'r')
    end
    plot(xlim,[true_rho true_rho])
    title(sprintf('%1.3f',pval))
    pause
    close(gcf)
end


end