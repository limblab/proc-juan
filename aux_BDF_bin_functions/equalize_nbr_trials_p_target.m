%
% Equalize nbr trial repetitions for all the targets
%
%   function single_trial_data = equalize_nbr_trials_p_target( single_trial_data )
%
%
% Note:
%   - not updating a bunch of meta variables



function single_trial_data = equalize_nbr_trials_p_target( single_trial_data )


% get the number of trials per target (excluding the last one, which is all
% of the previous one concatenated)

nbr_trials_p_target     = cellfun( @(x) size(x.neural_data.fr,3), single_trial_data.target );

nbr_trials_p_target     = nbr_trials_p_target(1:end-1);
min_nbr_trials_p_tgt    = min(nbr_trials_p_target);

nbr_targets             = length(single_trial_data.target)-1;

% -------------------------------------------------------------------------
% get rid of the extra trials
for t = 1:nbr_targets
    
    % neural data
    single_trial_data.target{t}.neural_data.fr = ...
        single_trial_data.target{t}.neural_data.fr(:,:,1:min_nbr_trials_p_tgt);
    single_trial_data.target{t}.neural_data.smoothed_fr = ...
        single_trial_data.target{t}.neural_data.smoothed_fr(:,:,1:min_nbr_trials_p_tgt);
    
    % dim_red neural data
    if isfield(single_trial_data.target{t}.neural_data,'dim_red')
        single_trial_data.target{t}.neural_data.dim_red.st_scores = ...
            single_trial_data.target{t}.neural_data.dim_red.st_scores(:,:,1:min_nbr_trials_p_tgt);
    end
    
    
    % emg data
    single_trial_data.target{t}.emg_data.emg = ...
        single_trial_data.target{t}.emg_data.emg(:,:,1:min_nbr_trials_p_tgt);
    
    % dim red emg data
    if isfield(single_trial_data.target{t}.emg_data,'dim_red')
        single_trial_data.target{t}.emg_data.dim_red.st_scores = ...
            single_trial_data.target{t}.emg_data.dim_red.st_scores(:,:,1:min_nbr_trials_p_tgt);
    end
    
    
    % pos data
    if isfield(single_trial_data.target{t},'pos')
        single_trial_data.target{t}.pos.data = ...
            single_trial_data.target{t}.pos.data(:,:,1:min_nbr_trials_p_tgt);
    end
    
    % vel data
    if isfield(single_trial_data.target{t},'vel')
        single_trial_data.target{t}.vel.data = ...
            single_trial_data.target{t}.vel.data(:,:,1:min_nbr_trials_p_tgt);
    end
    
    % force data
    if isfield(single_trial_data.target{t},'force')
        warning('force data not implemented yet');
    end
end


% -------------------------------------------------------------------------
% now fill the last trial, which is all the concatenated trials

aux_fr                  = single_trial_data.target{1}.neural_data.fr;
aux_fr_m                = single_trial_data.target{1}.neural_data.mn;
aux_fr_sd               = single_trial_data.target{1}.neural_data.sd;

if isfield(single_trial_data.target{1}.neural_data,'smoothedspikerate')
    aux_smoothed_fr     = single_trial_data.target{1}.neural_data.smoothed_fr;
    aux_smoothed_fr_m   = single_trial_data.target{1}.neural_data.smoothed_fr_mn;
    aux_smoothed_fr_sd  = single_trial_data.target{1}.neural_data.smoothed_fr_sd;
end

aux_emg                 = single_trial_data.target{1}.emg_data.emg;
aux_emg_m               = single_trial_data.target{1}.emg_data.mn;
aux_emg_sd              = single_trial_data.target{1}.emg_data.sd;

if isfield(single_trial_data.target{1},'pos')
    aux_pos             = single_trial_data.target{1}.pos.data;
    aux_pos_m           = single_trial_data.target{1}.pos.mn;
    aux_pos_sd          = single_trial_data.target{1}.pos.sd;
    aux_vel             = single_trial_data.target{1}.vel.data;
    aux_vel_m           = single_trial_data.target{1}.vel.mn;
    aux_vel_sd          = single_trial_data.target{1}.vel.sd;

end

if isfield(single_trial_data.target{1}.neural_data,'dim_red')
    aux_neural_scores   = single_trial_data.target{1}.neural_data.dim_red.st_scores;
    aux_neural_scores_m = single_trial_data.target{1}.neural_data.dim_red.st_scores_mn;
    aux_neural_scores_sd = single_trial_data.target{1}.neural_data.dim_red.st_scores_sd;
end

if isfield(single_trial_data.target{1}.emg_data,'dim_red')
    aux_emg_scores      = single_trial_data.target{1}.emg_data.dim_red.st_scores;
    aux_emg_scores_m    = single_trial_data.target{1}.emg_data.dim_red.st_scores_mn;
    aux_emg_scores_sd   = single_trial_data.target{1}.emg_data.dim_red.st_scores_sd;
end

% add the rest of the targets
for i = 2:nbr_targets
    
    aux_fr              = cat(3,aux_fr,single_trial_data.target{i}.neural_data.fr);
    aux_fr_m            = cat(1,aux_fr_m,single_trial_data.target{i}.neural_data.mn);
    aux_fr_sd           = cat(1,aux_fr_sd,single_trial_data.target{i}.neural_data.sd);
    
    if isfield(single_trial_data.target{1}.neural_data,'smoothedspikerate') 
        aux_smoothed_fr = cat(3,aux_smoothed_fr,single_trial_data.target{i}.neural_data.smoothed_fr);
        aux_smoothed_fr_m = cat(3,aux_smoothed_fr,single_trial_data.target{i}.neural_data.smoothed_fr_mn);
        aux_smoothed_fr_sd = cat(3,aux_smoothed_fr,single_trial_data.target{i}.neural_data.smoothed_fr_sd);
    end
    
    aux_emg             = cat(3,aux_emg,single_trial_data.target{i}.emg_data.emg);
    aux_emg_m           = cat(1,aux_emg_m,single_trial_data.target{i}.emg_data.mn);
    aux_emg_sd          = cat(1,aux_emg_sd,single_trial_data.target{i}.emg_data.sd);
    
    if isfield(single_trial_data.target{1},'pos')
        aux_pos         = cat(3,aux_pos,single_trial_data.target{i}.pos.data);
        aux_pos_m       = cat(1,aux_pos_m,single_trial_data.target{i}.pos.mn);
        aux_pos_sd      = cat(1,aux_pos_sd,single_trial_data.target{i}.pos.sd);
        aux_vel         = cat(3,aux_vel,single_trial_data.target{i}.vel.data);
        aux_vel_m       = cat(1,aux_vel_m,single_trial_data.target{i}.vel.mn);
        aux_vel_sd      = cat(1,aux_vel_sd,single_trial_data.target{i}.vel.sd);
    end
    
    if isfield(single_trial_data.target{1}.neural_data,'dim_red')
        aux_neural_scores   = cat(3,aux_neural_scores,single_trial_data.target{i}.neural_data.dim_red.st_scores);
        aux_neural_scores_m = cat(1,aux_neural_scores_m,single_trial_data.target{i}.neural_data.dim_red.st_scores_mn);
        aux_neural_scores_sd = cat(1,aux_neural_scores_sd,single_trial_data.target{i}.neural_data.dim_red.st_scores_sd);
    end
    
    if isfield(single_trial_data.target{1}.emg_data,'dim_red')
        aux_emg_scores      = cat(3,aux_emg_scores,single_trial_data.target{i}.emg_data.dim_red.st_scores);
        aux_emg_scores_m    = cat(1,aux_emg_scores_m,single_trial_data.target{i}.emg_data.dim_red.st_scores_mn);
        aux_emg_scores_sd   = cat(1,aux_emg_scores_sd,single_trial_data.target{i}.emg_data.dim_red.st_scores_sd);
    end
end


% now replace the last target in the struct by these vars

single_trial_data.target{end}.neural_data.fr     = aux_fr;
single_trial_data.target{end}.neural_data.mn     = aux_fr_m;
single_trial_data.target{end}.neural_data.sd     = aux_fr_sd;

if isfield(single_trial_data.target{1}.neural_data,'smoothedspikerate') 
    single_trial_data.target{end}.neural_data.smoothed_fr    = aux_smoothed_fr;
    single_trial_data.target{end}.neural_data.smoothed_fr_mn = aux_smoothed_fr_m;
    single_trial_data.target{end}.neural_data.smoothed_fr_sd = aux_smoothed_fr_sd;    
end

single_trial_data.target{end}.emg_data.emg       = aux_emg;
single_trial_data.target{end}.emg_data.mn        = aux_emg_m;
single_trial_data.target{end}.emg_data.sd        = aux_emg_sd;

if isfield(single_trial_data.target{1},'pos')
    single_trial_data.target{end}.pos.data       = aux_pos;
    single_trial_data.target{end}.pos.mn         = aux_pos_m;
    single_trial_data.target{end}.pos.sd         = aux_pos_sd;
    single_trial_data.target{end}.vel.data       = aux_vel;
    single_trial_data.target{end}.vel.mn         = aux_vel_m;
    single_trial_data.target{end}.vel.sd         = aux_vel_sd;
end

if isfield(single_trial_data.target{1}.neural_data,'dim_red')
    single_trial_data.target{end}.neural_data.dim_red.st_scores  = aux_neural_scores;
    single_trial_data.target{end}.neural_data.dim_red.mn         = aux_neural_scores_m;
    single_trial_data.target{end}.neural_data.dim_red.sd         = aux_neural_scores_sd;
end

if isfield(single_trial_data.target{1}.emg_data,'dim_red')
    single_trial_data.target{end}.emg_data.dim_red.st_scores     = aux_emg_scores;
    single_trial_data.target{end}.emg_data.dim_red.mn            = aux_emg_scores_m;
    single_trial_data.target{end}.emg_data.dim_red.sd            = aux_emg_scores_sd;
end