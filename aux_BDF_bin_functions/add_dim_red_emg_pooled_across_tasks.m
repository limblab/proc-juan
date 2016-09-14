%
% Computes muscle synergies with either PCA or NMF, and adds them to a
% struct of type datasets, both as a 'dim_red_emg' field, and within the
% single_trial_data substruct
%
%
% Inputs:
%   data_set            : struct of type datasets (as defined in
%                           batch_preprocess_dim_red_data)
%   params              : params for dim red; struct of type
%                           batch_preprocess_dim_red_data_params
%
% Outputs:
%   data_set            : struct of type datasets with new fields
%                           dim_red_emg, and dim_red_emg data in 
%
%

function data_set = add_dim_red_emg_pooled_across_tasks( data_set, params )


% -------------------------------------------------------------------------
% convert input struct params to input variables for
% dim_reduction_muscles_pooling_tasks(), which is the function that will do
% PCA/NMF of the EMGs

method              = params.dim_red_emg;
chosen_emgs         = data_set.chosen_emgs; % can pass it empty, as the single trial data that will be used already only has this nbr of channels
labels              = data_set.labels;
nbr_factors         = params.emg_factors;
plot_yn             = false;


% do dimensionality reduction
dim_red_emg         = dim_reduction_muscles_pooling_tasks( data_set, method, ...
                            chosen_emgs, nbr_factors, plot_yn );


% -------------------------------------------------------------------------
% now store the data in: 1) data_set.dim_red_emg; 2)
% data_set.stdata.emg_data.dim_red

% -------------
% 1. create dim_red_emg field in data_set
data_set.dim_red_emg = dim_red_emg;

% -------------
% 2. create a dim_red field for every target in
% data_set.stdata{i}.target{jj}.emg_data
%
% Note that each target in each stdata struct has the positions of the bins
% that correspond to each trial, in target{t}.bin_indx_p_trial (and also
% concatenated in bin_indx). The total number of bins per task is also
% known, so that can be used to save the dim_red emg data per trial


nbr_bdfs            = length(data_set.stdata);

% get the total nbr of bins per task, which will be used to store the dim
% red emgs (i.e. the scores)
nbr_bins_p_task     = cell2mat(cellfun( @(x) numel(x.target{end}.bin_indx), ...
                        data_set.stdata, 'UniformOutput', false ));

% get nbr targets per task
nbr_targets_p_task  = cell2mat(cellfun( @(x) numel(x.target)-1, data_set.stdata, ...
                        'UniformOutput', false)); 

% get nbr bins per target, for each task
nbr_bins_p_target_all_tasks = zeros(nbr_bdfs,max(nbr_targets_p_task));
for i = 1:nbr_bdfs
    aux_nbr_bins    = cellfun( @(x) length(x.conc_t), data_set.stdata{i}.target );
    nbr_bins_p_target_all_tasks(i,1:length(aux_nbr_bins)-1) = aux_nbr_bins(1:end-1);
    clear aux_nbr_bins
end


% create a bins_offset variable, to know by how many bins one has to offset
% the indexes depending on the task number
bins_offset_p_task  = [0, nbr_bins_p_task(1:end-1)];
                    
% the same, depending on the target nbr
bins_offset_p_target_n_task = [zeros(nbr_bdfs,1), nbr_bins_p_target_all_tasks(:,1:end-1)];

% check that the total nbr of bins per task equals the number of bins in
% dim_red_emg, 
if size(dim_red_emg.scores,1) ~= sum(nbr_bins_p_task)
    error('pooled_dim_red_emg: something went wrong, the total number of bins does not match the single trial data');
end


% do for each task 
for f = 1:length(data_set.stdata)
    
    % for each individual target 
    for t = 1:nbr_targets_p_task(f)
        
        
        % -----------------------------------------------------------------
        % ------------------------------
        % preallocate the vars
        data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores = ...
            zeros( size(data_set.stdata{f}.target{t}.bin_indx_p_trial,1), ...
                nbr_factors, size(data_set.stdata{f}.target{t}.bin_indx_p_trial,2) );
            
        data_set.stdata{f}.target{t}.emg_data.dim_red.scores = ...
            zeros( size(data_set.stdata{f}.target{t}.emg_data.conc_emg,1), ...
                nbr_factors );
        
        stat_fields = {'st_scores_mn','st_scores_sd'};
        for s = 1:length(stat_fields)
            data_set.stdata{f}.target{t}.emg_data.dim_red.(stat_fields{s}) = ...
                zeros( size(data_set.stdata{f}.target{t}.bin_indx_p_trial,1), ...
                    nbr_factors );
        end
          
        
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % save the data
        %
        % Note that to cut the single trials out of the scores we can't
        % trust the bin_indxs because those are relative to the raw
        % "cropped binned data", not to this ordered data
        
      
        % -----------------------------------------------------------------
        % a) the ("continuous") scores 
        aux_indx = bins_offset_p_task(f)+bins_offset_p_target_n_task(f,t) + ...
                    [1:length(data_set.stdata{f}.target{t}.conc_t)];        
        data_set.stdata{f}.target{t}.emg_data.dim_red.scores = ...
            dim_red_emg.scores(aux_indx,:);
        
        
        % -----------------------------------------------------------------
        % b) the scores for each trial
        nbr_trials_this     = size(data_set.stdata{f}.target{t}.bin_indx_p_trial,2);
        aux_indx            = reshape(aux_indx,[],nbr_trials_this);
        for r = 1:size(aux_indx,2)
            data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores(:,:,r) = ...
                dim_red_emg.scores(aux_indx(:,r),:);
        end
        
        % -----------------------------------------------------------------
        % c) the mean and SD
        data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores_mn = ...
            mean(data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores,3);
        data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores_sd = ...
            std(data_set.stdata{f}.target{t}.emg_data.dim_red.st_scores,0,3);
    end
    
    
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------   
    % for the concatenated (last) target
        
    % 1) some info about PCA/NMF, if it is the last target
    data_set.stdata{f}.target{end}.emg_data.dim_red.method = ...
        params.dim_red_emg;

    data_set.stdata{f}.target{end}.emg_data.dim_red.w = ...
        dim_red_emg.w;

    
    % -----------------------------------------------------------------
    % 2) the channels, for consistency
    data_set.stdata{f}.target{end}.emg_data.dim_red.chs = ...
        chosen_emgs;

    % -----------------------------------------------------------------
    % 3) the scores, and the mean / SD
    %
    % for the last target (which is all the concatenated targets) the
    % mean and SD are the concatenated mean SD for all the trials
            
    % vars to update
    these_vars = {'scores','st_scores_mn','st_scores_sd'};
    for v = 1:length(these_vars)
        % store the data from target 1
        data_set.stdata{f}.target{end}.emg_data.dim_red.(these_vars{v}) = ...
            data_set.stdata{f}.target{1}.emg_data.dim_red.(these_vars{v});
        for e = 2:nbr_targets_p_task(f)
            data_set.stdata{f}.target{end}.emg_data.dim_red.(these_vars{v}) = ...
                cat(1,data_set.stdata{f}.target{end}.emg_data.dim_red.(these_vars{v}),...
                data_set.stdata{f}.target{e}.emg_data.dim_red.(these_vars{v}));
        end
    end
    
    % -----------------------------------------------------------------
    % 3) the scores for each trial
    
    data_set.stdata{f}.target{end}.emg_data.dim_red.st_scores = ...
        data_set.stdata{f}.target{1}.emg_data.dim_red.st_scores;
    for e = 2:nbr_targets_p_task(f)
        data_set.stdata{f}.target{end}.emg_data.dim_red.st_scores = ...
            cat(3,data_set.stdata{f}.target{end}.emg_data.dim_red.st_scores,...
            data_set.stdata{f}.target{e}.emg_data.dim_red.st_scores);
    end
end

