%
% combine sorted neurons from the same channel of a TrialData file
%

function trial_data = mergeSortedNeurons( trial_data )

% check what arrays are included in the data
unit_guides         = getTDfields(trial_data,'unit_guides');
spike_names         = getTDfields(trial_data,'spikes');

% find the units that correspond to each channel in each array
for a = 1:length(spike_names)
    % find the channels in each array
    aux_units       = trial_data.(unit_guides{a});
    chs             = unique(aux_units(:,1));
    % save where the channels are in the unit matrices
    for c = 1:length(chs)
        indx_chs{a}.ch{c} = find( aux_units(:,1) == chs(c) );
    end
end

% do the merging
for i = 1:length(trial_data)
    for a = 1:length(spike_names)
        % fill aux_elecs with the number of spikes per channel in each bin
        % of this trial
        aux_elecs   = zeros(size(trial_data(i).(spike_names{a}),1),length(indx_chs{a}.ch));
        for c = 1:length(indx_chs{a}.ch)
            aux_elecs(:,c) = sum( trial_data(i).(spike_names{a})(:,indx_chs{a}.ch{c}), 2 );
        end    
        % overwrite the ARRAY_spikes matrix with sorted neurons with the
        % unsorted
        trial_data(i).(spike_names{a}) = aux_elecs;
    end
end