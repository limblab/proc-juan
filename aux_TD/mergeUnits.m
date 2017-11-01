%
% combine sorted neurons from the same channel of a TrialData file
%

function trial_data = mergeUnits( trial_data )


% figure out which arrays are here
arrays              = getTDfields(trial_data,'arrays');


% check what arrays are included in the data
unit_guides         = getTDfields(trial_data,'unit_guides');
spike_names         = getTDfields(trial_data,'spikes');


for idx_array = 1:length(arrays)
   
    array           = arrays{idx_array};
    
    master_sg       = trial_data(1).([array '_unit_guide']);

    % get rid of electrodes >96 (we only use UEAs with this many channels
    master_sg(master_sg(:,1)>96,:) = [];
    
    % get what electrodes we have
    elecs           = unique(master_sg(:,1));
    
    % and create new unit guide to overwrite it
    new_master_sg   = [elecs, zeros(length(elecs),1)];
    
    % Do: Merge units for each trial
    for i = 1:length(trial_data)
        
        temp_data   = trial_data(i).([array '_spikes']);
        
        merged_data = zeros(size(trial_data(i).([array '_spikes']),1),length(elecs));

        % merge each electrode
        for e = length(elecs)
            temp_spikes = temp_data(:,master_sg(:,1)==elecs(e));
            sum_temp = sum(temp_spikes,2);
            
            merged_data(:,e) = sum_temp;
        end
        
        % save unsorted units and update unit guide
        trial_data(i).([array '_spikes']) = merged_data;
        trial_data(i).([array '_unit_guide']) = new_master_sg;
    end
end