%
% Add dimensionally reduced neural and EMG data to a single_trial_data
% struct.
%
%   function single_trial_data = add_dim_red_data_to_single_trial_data( ...
%                               single_trial_data, varargin )
%
% Inputs (opt)          : [default]
%   single_trial_data   : single_trial_data struct, to which the data will
%                           be added
%   neural_data         : [empty] the dimensionally-reduced neural data
%                           that will be added. It has to be a dim_red_FR
%                           struct (generated with dim_reduction.m) 
%   emg_data            : [empty] the dimensionally-reduced emg data that
%                           will be added. It has to be a dim_red_emg
%                           struct (generated with dim_reduction_muscles.m)
% 
% Outputs:
%   single_trial_data   : the input struct with the added fields
%
%
%
% Note: 
%   - Not adding the "real time" to the time in dim_red_emg /
%   dim_red_neural structs, because it is measured w.r.t
%   cropped_binned_data (rather than binned_data), which means that targets
%   may start back to back  
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
    nbr_bins            = size(single_trial_data.target{1}.bin_indx_p_trial,1);
    nbr_chs             = length(neural_data.eigen);
    
    % -----------------------------------
    % 1) add scores (projections onto PCs) for each target
    
    for t = 1:nbr_targets
        
        % individual trials
        nbr_trials_this_tgt = size(single_trial_data.target{t}.bin_indx_p_trial,2);
        
        single_trial_data.target{t}.neural_data.dim_red.st_scores = ...
            zeros(nbr_bins,nbr_chs,nbr_trials_this_tgt);

        for r = 1:nbr_trials_this_tgt
            single_trial_data.target{t}.neural_data.dim_red.st_scores(:,:,r) = ...
                neural_data.scores(single_trial_data.target{t}.bin_indx_p_trial(:,r),:);
        end
            
        % all concatenated trials
        single_trial_data.target{t}.neural_data.dim_red.scores = ...
            neural_data.scores(single_trial_data.target{t}.bin_indx,:); 
%         single_trial_data.target{t}.neural_data.dim_red.t = ...
%             neural_data.t(single_trial_data.target{t}.bin_indx);
        
        % add mean and SD scores
        single_trial_data.target{t}.neural_data.dim_red.st_scores_mn = ...
            mean(single_trial_data.target{t}.neural_data.dim_red.st_scores,3);
        single_trial_data.target{t}.neural_data.dim_red.st_scores_sd = ...
            std(single_trial_data.target{t}.neural_data.dim_red.st_scores,0,3);
    end
    
    % -----------------------------------
    % 2) fill last field of single_trial_data, which are all the
    % concatenated trials
    
    % init vectors to concatenate data
    aux_st_scores       = single_trial_data.target{1}.neural_data.dim_red.st_scores;
    aux_scores          = single_trial_data.target{1}.neural_data.dim_red.scores;
%     aux_t               = single_trial_data.target{1}.neural_data.dim_red.t;
    
    aux_st_scores_mn    = single_trial_data.target{1}.neural_data.dim_red.st_scores_mn;
    aux_st_scores_sd    = single_trial_data.target{1}.neural_data.dim_red.st_scores_sd;
    
    % add info for all the concatenated targets, as well as the scores
    single_trial_data.target{end}.neural_data.dim_red.method = ...
        neural_data.method;
    single_trial_data.target{end}.neural_data.dim_red.w = ...
        neural_data.w;
    
    % concatenate scores & time
    for t = 2:nbr_targets
        aux_st_scores  	= cat(3,aux_st_scores,single_trial_data.target{t}.neural_data.dim_red.st_scores);
        aux_scores      = cat(1,aux_scores,single_trial_data.target{t}.neural_data.dim_red.scores);
%         aux_t           = cat(1,aux_t,single_trial_data.target{t}.neural_data.dim_red.t);
        
        aux_st_scores_mn = cat(1,aux_st_scores_mn,single_trial_data.target{t}.neural_data.dim_red.st_scores_mn);
        aux_st_scores_sd = cat(1,aux_st_scores_sd,single_trial_data.target{t}.neural_data.dim_red.st_scores_sd);
    end
%     single_trial_data.target{end}.neural_data.dim_red.t = aux_t;
    single_trial_data.target{end}.neural_data.dim_red.scores = aux_scores;
    
    single_trial_data.target{end}.neural_data.dim_red.eigen = ...
        neural_data.eigen;
    single_trial_data.target{end}.neural_data.dim_red.chs = ...
        neural_data.chs;    
    
    single_trial_data.target{end}.neural_data.dim_red.st_scores = ...
        aux_st_scores;
    single_trial_data.target{end}.neural_data.dim_red.st_scores_mn = ...
        aux_st_scores_mn;
    single_trial_data.target{end}.neural_data.dim_red.st_scores_sd = ...
        aux_st_scores_sd;
end


% -------------------------------------------------------------------------
% add emg data
if emg_data_yn
    error ('adding dimensionally-reduced EMG not implemented yet')
end
