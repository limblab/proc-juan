%
% Wrapper to do dPCA. It is written to do dimensionality reduction
% comparing different tasks, so the covariates are hardcoded to be 'task',
% 'target', 'time', and 'task/target interaction'. Also, for this study we
% compared tasks with 6 or 8 targets so the code is also written to choose
% the targets in the 8 target task that best match the targets in the
% 6-target tasks.
%

function dPCA_results = call_dPCA( single_trial_data, varargin )


% read input params
if nargin >= 2
    num_comps       = varargin{1};
else
    num_comps       = 15;
end
if nargin == 3
    plot_yn         = varargin{2};
else
    plot_yn         = true;
end
    
% retrieve neural chs
neural_chs          = single_trial_data{1}.target{1}.neural_data.neural_chs;

    
% check if the dimensions in single_trial_data are consistent
if size(single_trial_data{1}.target{1}.neural_data.smoothed_fr,2) ~= numel(neural_chs)
    error('single trial data does not include all the neural channels')
end

% make the target averaged responses for each task have equal length.
% This is not done in single_trial_analysis.m, where single trial
% duration is only equalized for each task
single_trial_data   = equalize_single_trial_dur( single_trial_data );



% get rid of the last target, which is all the concatenated targets
for i = 1:length(single_trial_data)
    single_trial_data{i}.target(end) = [];
end

% ------------------------------------------------------------------------
% 1. arrange the data (as described in dpca_demo)

% N is the number of neurons
% S is the number of conditions --> tasks in our case
% D is the number of decisions --> targets in our case
%       ToDo: so far we are choosing the min, but see if they can be different for each task
% T is the number of time points --each trial should have the same duration
% in time !!!
N                   = numel(neural_chs);
S                   = numel(single_trial_data);
D                   = min(cellfun(@(x) length(x.target), single_trial_data));
T                   = size(single_trial_data{1}.target{1}.neural_data.fr,1);
% max number of repetitions
max_trial_num       = 0;
for i = 1:S
    if max(cellfun(@(x) size(x.neural_data.fr,3), single_trial_data{1}.target )) > max_trial_num
        max_trial_num = max(cellfun(@(x) size(x.neural_data.fr,3), single_trial_data{1}.target ));
    end
end


% trial_num: N x S x D
trial_num           = zeros(N,S,D);

% firing_rates: N x S x D x T x max_trial_num -- populated with our
% single_trial_data
firing_rates        = nan(N,S,D,T,max_trial_num);
% ·········································································
% In the 1D tasks, targets 1 to 6 go from left to right, but in the 2D
% tasks targets are ordered as follows: 5, 7, 8, 6, 4, 2, 1, 3 --beginning
% at 12 o'clock and going clockwise. They will be paired as 1D/2D: 1/1,
% 2/2, 3/3, 4/6, 5/7, 6, 8 (as defined in target_order)
% ·········································································
target_order        = [1 2 3 4 5 6; 1 2 3 6 7 8]; % row 1: 1D task; row 2: 2D task

for n = 1:N
    for s = 1:S
        for d = 1:D
            if length(single_trial_data{s}.target) == D
                tgt = target_order(1,d);
            elseif length(single_trial_data{s}.target) >= D
                tgt = target_order(2,d);
            end
            trials_this = size(single_trial_data{s}.target{tgt}.neural_data.smoothed_fr(:,n,:),3);
            firing_rates(n,s,d,:,1:trials_this) = ...
                squeeze(single_trial_data{s}.target{tgt}.neural_data.smoothed_fr(:,n,:));
        end
    end
end


% firing_rates_average: N x S x D x T -- these are PSTHs
firing_rates_avg    = nanmean(firing_rates, 5);

% ------------------------------------------------------------------------
% 2. Define parameters

% 1 - task
% 2 - target
% 3 - time
% [1 3] - task/time interaction
% [2 3] - target/time interaction
% [1 2] - task/target interaction
% [1 2 3] - rest

% we have datasets that only have one target (D=1), like the ball task, in
% that case the params have to be different
if D > 1
    combined_params = { {1,[1 3]}, {2,[2,3]}, {3}, {[1 2],[1 2 3]} };
    marg_names      = {'task','target','time','task/target interaction'};
    marg_colors     = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256; % blue, red, grey, purple

    time            = (1:T)*single_trial_data{1}.target{1}.bin_size;
    time_events     = 1;
else
    combined_params = { {1,[1 3]}, {3} };
    marg_names      = {'task','time'};
    marg_colors     = [23 100 171; 187 20 25; 150 150 150]/256;

    time            = (1:T)*single_trial_data{1}.target{1}.bin_size;
    time_events     = 1;    
end
    

% ------------------------------------------------------------------------
% 2. Do dPCA without regularization


[W, V, which_marg]  = dpca( firing_rates_avg, num_comps, 'combinedParams', combined_params );

expl_var            = dpca_explainedVariance(firing_rates_avg, W, V, ...
                        'combinedParams', combined_params);

if plot_yn                     
    dpca_plot(firing_rates_avg, W, V, @dpca_plot_default_j, ...
        'explainedVar', expl_var, ...
        'marginalizationNames', marg_names, ...
        'marginalizationColours', marg_colors, ...
        'whichMarg', which_marg,                 ...
        'time', time,                        ...
        'timeEvents', time_events,               ...
        'timeMarginalization', 3, ...
        'legendSubplot', num_comps);                
end


% ------------------------------------------------------------------------
% 2. Do dPCA in each marginalization separately 

% dpca_perMarginalization(firing_rates_avg, @dpca_plot_default, ...
%    'combinedParams', combined_params);


% ------------------------------------------------------------------------
% project data onto dPC axes
[lat_vars, lat_vars_st] = get_lat_vars_dPCA( firing_rates_avg, firing_rates, W, ...
    'explainedVar', expl_var, ...
    'marginalizationNames', marg_names, ...
    'marginalizationColours', marg_colors, ...
    'whichMarg', which_marg,                 ...
    'time', time,                        ...
    'timeEvents', time_events,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', num_comps);


% ------------------------------------------------------------------------
% 3. Return dPCA results

dPCA_results.W          = W;
dPCA_results.V          = V;
dPCA_results.lat_vars_mn = lat_vars;
dPCA_results.lat_vars   = lat_vars_st;
dPCA_results.which_marg = which_marg;
dPCA_results.marg_names = marg_names;
dPCA_results.expl_var   = expl_var;
dPCA_results.num_comps  = num_comps;
dPCA_results.combined_params = combined_params;