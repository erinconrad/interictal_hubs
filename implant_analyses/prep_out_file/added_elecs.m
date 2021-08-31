function [new_elecs,new_depths,really_unchanged] = added_elecs(name,chLabels,tentative_unchanged)

switch name
    case 'HUP099' 
        elecs = {'DHI','DHT','DHAA','DA'};
        depths = {'DHI','DHT','DHAA','DA'};
    case 'HUP100'
        elecs = {'LSFP','LIFP','LPOF','LOC','DSC','DSEP','DMPF','DAMF'};
        depths = {'DSC','DSEP','DMPF','DAMF'};
    case 'HUP110'
        elecs = {'LMF','LPF','LG'};
        depths = {};
    case 'HUP111'
        elecs = {'RG','ROF','RAST','RMST','RPST','RSP','RO','LAST','LPST'};
        depths = {};
    case 'HUP128'
        elecs = {'LW','LX','LY','LZ'};
        depths=  {'LW','LX','LY','LZ'};
    case 'HUP129'
        elecs = {'RX','RY','RZ'};
        depths = {'RX','RY','RZ'};
    case 'HUP132'
        elecs = {'LK','LM','LN'};
        depths = {'LK','LM','LN'};
    case 'HUP136'
        elecs = {'LM','LT','LI','LO'};
        depths = {'LM','LT','LI','LO'};
    case 'HUP152'
        elecs = {'LH','LP','LQ'};
        depths = {'LH','LP','LQ'};
    case 'HUP168'
        elecs = {'LN','LO','LP','LR'};
        depths = {'LN','LO','LP','LR'};
    case 'HUP193'
        elecs = {'LK','LL'};
        depths = {'LK','LL'};
    case 'HUP201'
        elecs = {'LM'};
        depths = {'LM'};
    case 'HUP209'
        elecs = {'LM','LN'};
        depths = {'LM','LN'};
end
    
    
    
%% Match up electrodes with these
elecs_added = zeros(length(chLabels),1);
depths_added = zeros(length(chLabels),1);
really_unchanged = ones(length(tentative_unchanged),1);
for i = 1:length(chLabels)
    for j = 1:length(elecs)
        if contains(chLabels{i},elecs{j})
            elecs_added(i) = 1;
        end
        
    end
    
    for j = 1:length(depths)
        if contains(chLabels{i},depths{j})
            depths_added(i) = 1;
        end
        
    end
    
end

for i = 1:length(really_unchanged)
    for j = 1:length(elecs)
        if contains(tentative_unchanged{i},elecs{j})
            really_unchanged(i) = 0;
        end
        
    end
    
    
end

%% Get corresponding labels
new_elecs = chLabels(logical(elecs_added));
new_depths = chLabels(logical(depths_added));
really_unchanged = tentative_unchanged(logical(really_unchanged));