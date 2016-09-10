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

% nbr of vars to add (1/2)
vars_to_add             = [neural_data_yn, emg_data_yn];
indx_vars_to_add        = find(vars_to_add);

% -------------------------------------------------------------------------
% add data

for v = 1:length(indx_vars_to_add)
    
    % if it is the neural data that we want to add
    if indx_vars_to_add(v) == 1
       
        this_var        = 'neural_data';
        nbr_chs         = length(neural_data.eigen);
        
    % if it is the emg data that we want to add
    elseif indx_vars_to_add(v) == 2

        this_var        = 'emg_data';
        nbr_chs         = length(emg_data.chs);
    end
    
    % get some meta info
    nbr_targets         = length(single_trial_data.target) - 1;
    nbr_bins            = size(single_trial_data.target{1}.bin_indx_p_trial,1);
    

    % -----------------------------------
    % 1) add scores (projections onto PCs) for each target

    for t = 1:nbr_targets

        % individual trials
        nbr_trials_this_tgt = size(single_trial_data.target{t}.bin_indx_p_trial,2);

        switch this_var
            case 'neural_data'
                single_trial_data.target{t}.(this_var).dim_red.st_scores = ...
                    zeros(nbr_bins,nbr_chs,nbr_trials_this_tgt);
            case 'emg_data'
                single_trial_data.target{t}.(this_var).dim_red.st_scores = ...
                    zeros(nbr_bins,size(emg_data.w,2),nbr_trials_this_tgt);
        end
        
        for r = 1:nbr_trials_this_tgt
            switch this_var
                case 'neural_data'
                    single_trial_data.target{t}.(this_var).dim_red.st_scores(:,:,r) = ...
                        neural_data.scores(single_trial_data.target{t}.bin_indx_p_trial(:,r),:);
                case 'emg_data'
                    single_trial_data.target{t}.(this_var).dim_red.st_scores(:,:,r) = ...
                        emg_data.scores(single_trial_data.target{t}.bin_indx_p_trial(:,r),:);
            end
        end

        % all concatenated trials
        switch this_var
            case 'neural_data'
                single_trial_data.target{t}.(this_var).dim_red.scores = ...
                    neural_data.scores(single_trial_data.target{t}.bin_indx,:); 
%         single_trial_data.target{t}.(this_var).dim_red.t = ...
%             (this_var).t(single_trial_data.target{t}.bin_indx);
            case 'emg_data'
                single_trial_data.target{t}.(this_var).dim_red.scores = ...
                    emg_data.scores(single_trial_data.target{t}.bin_indx,:); 
        end

        % add mean and SD scores
        single_trial_data.target{t}.(this_var).dim_red.st_scores_mn = ...
            mean(single_trial_data.target{t}.(this_var).dim_red.st_scores,3);
        single_trial_data.target{t}.(this_var).dim_red.st_scores_sd = ...
            std(single_trial_data.target{t}.(this_var).dim_red.st_scores,0,3);
    end

    % -----------------------------------
    % 2) fill last field of single_trial_data, which are all the
    % concatenated trials

    % init vectors to concatenate data
    aux_st_scores       = single_trial_data.target{1}.(this_var).dim_red.st_scores;
    aux_scores          = single_trial_data.target{1}.(this_var).dim_red.scores;
%     aux_t               = single_trial_data.target{1}.(this_var).dim_red.t;

    aux_st_scores_mn    = single_trial_data.target{1}.(this_var).dim_red.st_scores_mn;
    aux_st_scores_sd    = single_trial_data.target{1}.(this_var).dim_red.st_scores_sd;

    % add info for all the concatenated targets, as well as the scores
    switch this_var
        case 'neural_data'
            single_trial_data.target{end}.(this_var).dim_red.method = ...
                neural_data.method;
            single_trial_data.target{end}.(this_var).dim_red.w = ...
                neural_data.w;
        case 'emg_data'
            single_trial_data.target{end}.(this_var).dim_red.method = ...
                emg_data.method;
            single_trial_data.target{end}.(this_var).dim_red.w = ...
                emg_data.w;
    end

    % concatenate scores & time
    for t = 2:nbr_targets
        aux_st_scores  	= cat(3,aux_st_scores,single_trial_data.target{t}.(this_var).dim_red.st_scores);
        aux_scores      = cat(1,aux_scores,single_trial_data.target{t}.(this_var).dim_red.scores);
%         aux_t           = cat(1,aux_t,single_trial_data.target{t}.(this_var).dim_red.t);

        aux_st_scores_mn = cat(1,aux_st_scores_mn,single_trial_data.target{t}.(this_var).dim_red.st_scores_mn);
        aux_st_scores_sd = cat(1,aux_st_scores_sd,single_trial_data.target{t}.(this_var).dim_red.st_scores_sd);
    end
%     single_trial_data.target{end}.(this_var).dim_red.t = aux_t;
    single_trial_data.target{end}.(this_var).dim_red.scores = aux_scores;

    switch this_var
        case 'neural_data'
            single_trial_data.target{end}.(this_var).dim_red.eigen = ...
                neural_data.eigen;
            single_trial_data.target{end}.(this_var).dim_red.chs = ...
                neural_data.chs;    
        case 'emg_data'
            if isfield(emg_data,'eigen')
                single_trial_data.target{end}.(this_var).dim_red.eigen = ...
                    emg_data.eigen;
            end
            single_trial_data.target{end}.(this_var).dim_red.chs = ...
                emg_data.chs;    
    end

    single_trial_data.target{end}.(this_var).dim_red.st_scores = ...
        aux_st_scores;
    single_trial_data.target{end}.(this_var).dim_red.st_scores_mn = ...
        aux_st_scores_mn;
    single_trial_data.target{end}.(this_var).dim_red.st_scores_sd = ...
        aux_st_scores_sd;
end
