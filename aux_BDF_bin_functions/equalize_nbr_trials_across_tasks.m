%
% Equalize the number of trials across different tasks.
%
%   function single_trial_data = equalize_nbr_trials_across_tasks( single_trial_data, varargin )
%
%
% Inputs (opt)          : [default]
%   single_trial_data   : cell array of single_trial_data structs (one
%                           field per task)
%   (target)           	: ['all_conc'] 'all_conc': will cut the
%                           concatenated target (last field in struct) to
%                           the same length, also cutting the trials of
%                           each target to the corresponding number. ['X']
%                           (X=target number): will cut each of the trials
%                           to the same number of repetitions, also
%                           overwriting the last field.
% 
% Outputs:
%   single_trial_data   : single_trial_data struct with updated nbr of
%                           trials
%
%

function single_trial_data = equalize_nbr_trials_across_tasks( single_trial_data, varargin )


% -------------------------------------------------------------------------
% Read input params
if nargin == 2
    target              = varargin{1};
else
    target              = 'all_conc';
end


% -------------------------------------------------------------------------
% get some info


% number of tasks
nbr_bdfs                = length(single_trial_data);

% get the number of targets per task
nbr_targets_p_task      = cellfun( @(x) length(x.target)-1, single_trial_data );

% get the number of trials per target
nbr_trials_p_target_n_task  = cellfun( @(x) size(x.target{1}.neural_data.fr,3), ...
                                single_trial_data );

% and total number of targets of that task
total_nbr_trials_p_task = cellfun( @(x) size(x.target{end}.neural_data.fr,3), ...
                                single_trial_data );


% -------------------------------------------------------------------------
% cut to the selected number of trials, depending on what we want to
% compare

switch target
    case 'all_conc'
        target_to_edit  = nbr_targets_p_task + 1;
        nbr_trials_to_keep  = min(total_nbr_trials_p_task);
        % get number of trials to keep
        nbr_trials_p_tgt_to_keep = repmat(nbr_trials_to_keep,1,2)./nbr_targets_p_task;
    otherwise
        target_to_edit  = target;
        nbr_trials_p_tgt_to_keep = min(nbr_trials_p_target_n_task);
        warning('this mode has not been tested');
        pause;
end


% -------------------------------------------------------------------------
% do for each task


% get field names, for updating the number of targets
neural_names            = fieldnames(single_trial_data{1}.target{end}.neural_data);
emg_names               = fieldnames(single_trial_data{1}.target{end}.emg_data);
if isfield(single_trial_data{1}.target{1}.neural_data,'dim_red')
%    neural_dim_red_names = fieldnames(single_trial_data{1}.target{end}.neural_data.dim_red);
    neural_dim_red_names = fieldnames(single_trial_data{1}.target{1}.neural_data.dim_red);
end
if isfield(single_trial_data{1}.target{1}.emg_data,'dim_red')
    emg_dim_red_names   = fieldnames(single_trial_data{1}.target{end}.emg_data.dim_red);
end
if isfield(single_trial_data{1}.target{1},'pos')
    pos_names           = fieldnames(single_trial_data{1}.target{end}.pos);
end
if isfield(single_trial_data{1}.target{1},'vel')
    vel_names           = fieldnames(single_trial_data{1}.target{end}.vel);
end
if isfield(single_trial_data{1}.target{1},'force')
    warning('EQUALIZE_SINGLE_TRIAL_DUR: force data length equalization not implemented yet');
end


for i = 1:2

    % -----------------
    % neural data
    for ii = 1:numel(neural_names)
        % see if it is a 3D matrix, because dim 3 is the trial nbr
        if size(single_trial_data{i}.target{1}.neural_data.(neural_names{ii}),3) > 1
            % do for each of the targets
            for t = 1:nbr_targets_p_task(i)
                single_trial_data{i}.target{t}.neural_data.(neural_names{ii}) = ...
                    single_trial_data{i}.target{t}.neural_data.(neural_names{ii})...
                    (:,:,1:nbr_trials_p_tgt_to_keep);
            end
            % and for the concatenated target (last target in the struct)
            aux_cat     = single_trial_data{i}.target{1}.neural_data.(neural_names{ii});
            for t = 2:nbr_targets_p_task(i)
                aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.neural_data.(neural_names{ii}));
            end
            single_trial_data{i}.target{end}.neural_data.(neural_names{ii}) = aux_cat;
        end
    end
        
	% -----------------
    % EMG data
    for ii = 1:numel(emg_names)
        % see if it is a 3D matrix, because dim 3 is the trial nbr
        if size(single_trial_data{i}.target{1}.emg_data.(emg_names{ii}),3) > 1
            % do for each of the targets
            for t = 1:nbr_targets_p_task(i)
                single_trial_data{i}.target{t}.emg_data.(emg_names{ii}) = ...
                    single_trial_data{i}.target{t}.emg_data.(emg_names{ii})...
                    (:,:,1:nbr_trials_p_tgt_to_keep);
            end
            % and for the concatenated target (last target in the struct)
            aux_cat     = single_trial_data{i}.target{1}.emg_data.(emg_names{ii});
            for t = 2:nbr_targets_p_task(i)
                aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.emg_data.(emg_names{ii}));
            end
            single_trial_data{i}.target{end}.emg_data.(emg_names{ii}) = aux_cat;
        end
    end
                    
    % -------------------
    % dim reduced neural data
    if exist('neural_dim_red_names','var')
         for ii = 1:numel(neural_dim_red_names)
            % see if it is a 3D matrix, because dim 3 is the trial nbr
            if size(single_trial_data{i}.target{1}.neural_data.dim_red.(neural_dim_red_names{ii}),3) > 1
                % do for each of the targets
                for t = 1:nbr_targets_p_task(i)
                    single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{ii}) = ...
                        single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{ii})...
                        (:,:,1:nbr_trials_p_tgt_to_keep);
                end
                % and for the concatenated target (last target in the struct)
                aux_cat     = single_trial_data{i}.target{1}.neural_data.dim_red.(neural_dim_red_names{ii});
                for t = 2:nbr_targets_p_task(i)
                    aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{ii}));
                end
                single_trial_data{i}.target{end}.neural_data.dim_red.(neural_dim_red_names{ii}) = aux_cat;
            end
        end
    end
    
    % -------------------
    % dim reduced emg data
    if exist('emg_dim_red_names','var')
         for ii = 1:numel(emg_dim_red_names)
            % see if it is a 3D matrix, because dim 3 is the trial nbr
            if size(single_trial_data{i}.target{1}.emg_data.dim_red.(emg_dim_red_names{ii}),3) > 1
                % do for each of the targets
                for t = 1:nbr_targets_p_task(i)
                    single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{ii}) = ...
                        single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{ii})...
                        (:,:,1:nbr_trials_p_tgt_to_keep);
                end
                % and for the concatenated target (last target in the struct)
                aux_cat     = single_trial_data{i}.target{1}.emg_data.dim_red.(emg_dim_red_names{ii});
                for t = 2:nbr_targets_p_task(i)
                    aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{ii}));
                end
                single_trial_data{i}.target{end}.emg_data.dim_red.(emg_dim_red_names{ii}) = aux_cat;
            end
        end
    end
    
    % -------------------
    % Pos data
    if exist('pos_names','var')
        for ii = 1:numel(pos_names)
            % see if it is a 3D matrix, because dim 3 is the trial nbr
            if size(single_trial_data{i}.target{1}.pos.(pos_names{ii}),3) > 1
                % do for each of the targets
                for t = 1:nbr_targets_p_task(i)
                    single_trial_data{i}.target{t}.pos.(pos_names{ii}) = ...
                        single_trial_data{i}.target{t}.pos.(pos_names{ii})...
                        (:,:,1:nbr_trials_p_tgt_to_keep);
                end
                % and for the concatenated target (last target in the struct)
                aux_cat     = single_trial_data{i}.target{1}.pos.(pos_names{ii});
                for t = 2:nbr_targets_p_task(i)
                    aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.pos.(pos_names{ii}));
                end
                single_trial_data{i}.target{end}.pos.(pos_names{ii}) = aux_cat;
            end
        end
    end
    
    % -------------------
    % Vel data
    if exist('vel_names','var')
        for ii = 1:numel(vel_names)
            % see if it is a 3D matrix, because dim 3 is the trial nbr
            if size(single_trial_data{i}.target{1}.vel.(vel_names{ii}),3) > 1
                % do for each of the targets
                for t = 1:nbr_targets_p_task(i)
                    single_trial_data{i}.target{t}.vel.(vel_names{ii}) = ...
                        single_trial_data{i}.target{t}.vel.(vel_names{ii})...
                        (:,:,1:nbr_trials_p_tgt_to_keep);
                end
                % and for the concatenated target (last target in the struct)
                aux_cat     = single_trial_data{i}.target{1}.vel.(vel_names{ii});
                for t = 2:nbr_targets_p_task(i)
                    aux_cat = cat(3,aux_cat,single_trial_data{i}.target{t}.vel.(vel_names{ii}));
                end
                single_trial_data{i}.target{end}.vel.(vel_names{ii}) = aux_cat;
            end
        end
    end
    
    % -------------------        
    % Update other stuff and meta info
    
end


% -------------------------------------------------------------------------
% overwrite the concatenated fields

for i = 1:nbr_bdfs
   
    single_trial_data{i} = concatenate_single_trials( single_trial_data{i} );
end

end
