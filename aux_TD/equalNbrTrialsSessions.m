%
% Equalize number of trials across sessions and targets
% 

function [trial_data, min_nbr_trials] = equalNbrTrialsSessions( trial_data )


% get the sessions
sessions    = unique({trial_data.date});
% get the targets
targets     = unique(cell2mat({trial_data.target_direction}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% Find minimum number of trials and targets across all session x target
% combinations


% get the number of trials per target and session
trials_p_target = zeros(length(sessions),length(targets));

for s = 1:length(sessions)
    for t = 1:length(targets)
        trials_p_target(s,t) = numel(getTDidx(trial_data,'date',sessions{s},...
                                'target_direction',targets(t)));
    end
end

% find the minimum number of trials across all sessions and targets
min_nbr_trials = min(min(trials_p_target));

disp(['Keeping the first ' num2str(min_nbr_trials) ' trials per target and session']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% Create the new trial_data struct


% preallocate new_trial_data -- dirty trick: i'm repeating the first field
% of trial_data to define new_trial_data, so they have the same structure
new_trial_data = repmat(trial_data(1),1,...
                    min_nbr_trials*length(sessions)*length(targets));

% populate new_trial_data
for s = 1:length(sessions)
    for t = 1:length(targets)
        % get the trial number for this session and target
        t_idx       = getTDidx(trial_data,'date',sessions{s},...
                        'target_direction',targets(t));
        % keep only the first min_nbr_trials
        t_idx       = t_idx(1:min_nbr_trials);
        
        % add them to the new_trial_data struct
        start_idx   = length(targets)*min_nbr_trials*(s-1) + min_nbr_trials*(t-1) + 1;
        end_idx     = start_idx + min_nbr_trials - 1;
        
        new_trial_data(start_idx:end_idx) = trial_data(t_idx);
    end
end

% overwrite trial_data to return

trial_data      = new_trial_data;