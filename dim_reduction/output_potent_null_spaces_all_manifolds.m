%
% Do output potent/null analysis of all tasks in a session. 
%
%
%

function opn_spaces = output_potent_null_spaces_all_manifolds( single_trial_data, varargin )


% read input params, if passed
if nargin == 2
    params              = output_potent_null_spaces_defaults( varargin{1} );
else
    params              = output_potent_null_spaces_defaults;
end


nbr_bdfs                = length(single_trial_data);

% -------------------------------------------------------------------------
% prepare data for the analysis

% 1) equalize trial duration across all tasks
single_trial_data       = equalize_single_trial_dur( single_trial_data, ...
                            'time_win', params.time_win );

% 2) equalize number of trials for all targets of a given task
for i = 1:nbr_bdfs
    single_trial_data{i} = equalize_nbr_trials_p_target( single_trial_data{i} );
end

% 3) equalize number of trials across tasks
single_trial_data       = equalize_nbr_trials_across_tasks( single_trial_data, params.target );


% -------------------------------------------------------------------------
% do !

for i = 1:nbr_bdfs
    opn_spaces(i)       = output_potent_null_spaces( single_trial_data{i}, params );
end