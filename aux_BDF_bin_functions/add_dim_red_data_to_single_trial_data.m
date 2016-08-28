%
% add dimensionally reduced neural and EMG data to a single trial data
% struct 
%

function single_trial_data = add_dim_red_data_to_single_trial_data( ...
                                single_trial_data, varargin )


% read input parameters
if nargin == 1
    neural_data         = [];
    emg_data            = [];
end
if nargin == 2
    neural_data         = varargin{1};
    emg_data            = [];
end
if nargin == 3
    neural_data         = varargin{1};
    emg_data            = varargin{2};
end

if ~isempty(neural_data)
    neural_data_yn      = true;
else
    neural_data_yn      = false;
end
if ~isempty(emg_data)
    emg_data_yn         = true;
else
    emg_data_yn         = false;
end


% -------------------------------------------------------------------------
% add neural data
if neural_data_yn
    
    nbr_targets         = length(single_trial_data.target) - 1;
    
    % -----------------------------------
    % 1) add scores (projections onto PCs) for each target
    
    for t = 1:nbr_targets
        single_trial_data.target{t}.neural_data.dim_red.scores = ...
            neural_data.scores(single_trial_data.target{t}.bin_indx,:); 
        single_trial_data.target{t}.neural_data.dim_red.t = ...
            neural_data.t(single_trial_data.target{t}.bin_indx);
    end
    
    % -----------------------------------
    % 2) fill last field of single_trial_data, which are all the
    % concatenated trials
    
    % init vectors to concatenate data
    aux_scores          = single_trial_data.target{1}.neural_data.dim_red.scores;
    aux_t               = single_trial_data.target{1}.neural_data.dim_red.t;
    
    % add info for all the concatenated targets, as well as the scores
    single_trial_data.target{end}.neural_data.dim_red.method = ...
        neural_data.method;
    single_trial_data.target{end}.neural_data.dim_red.w = ...
        neural_data.w;
    
    % concatenate scores & time
    for t = 2:nbr_targets
        aux_scores      = cat(1,aux_scores,single_trial_data.target{t}.neural_data.dim_red.scores);
        aux_t           = cat(1,aux_t,single_trial_data.target{t}.neural_data.dim_red.t);
    end
    single_trial_data.target{end}.neural_data.dim_red.t = aux_t;
    single_trial_data.target{end}.neural_data.dim_red.scores = aux_scores;
    
    single_trial_data.target{end}.neural_data.dim_red.eigen = ...
        neural_data.eigen;
    single_trial_data.target{end}.neural_data.dim_red.chs = ...
        neural_data.chs;
end


% -------------------------------------------------------------------------
% add emg data
if emg_data_yn
    error ('adding dimensionally-reduced EMG not implemented yet')
end
