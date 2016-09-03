% 
% Take an array of target averaged responses for several tasks obtained
% with single_trial_analysis_dim_red and make them all have the same
% duration (number of bins). This can be done in different ways:
%
%   function single_trial_data   = equalize_single_trial_dur( single_trial_data, varargin )
%
% Inputs (opt)              : [default]
%   single_trial_data       : single_trial_data struct, generated with
%                               get_single_trial_data.m
%   (mode)                  : ['min_dur'] method for equalizing trial
%                               duration: 1) 'min_dur': time warp to the
%                               duration of the shortest trials; 2)
%                               'max_dur': time warp to the duration of the
%                               longest trials; 3) 'time_win': start and
%                               end times of the window [t_init, t_end]
%                               with respect to how the trials have been
%                               cut.
%   (time_window)           : [t_init,t_end]: if mode == 'period', defines
%                               from when to when we'll keep the data p
% 
% Outputs:
%   single_trial_data       : single_trial_data with equalized trial
%                               duration for all tasks (and targets)
%
%
% ToDo's:
%   - need to implement updating the fields with information in
%   single_trial_data.m 
%   - need to update the time vector of the time warping option
%   - need to add force data


function single_trial_data   = equalize_single_trial_dur( single_trial_data, varargin )


% -------------------------------------------------------------------------
% Read input params
if nargin >= 2
    mode                = varargin{1};
else
    mode                = 'min_dur';
end
if strcmp(mode,'time_win')
    if nargin ~= 3, 
        error('Need to pass the time window for cutting, when mode == time'); 
    else
        time_window = varargin{2};
        if size(time_window) ~= [1 2]
            error('time window needs to have size [1,2] ([t_init,t_end])');
        end
    end
end


% -------------------------------------------------------------------------
% Decide to what number of bins to resample all trial averaged data


% get number of samples of the average responses for each dataset
nbr_bins_dataset        = cellfun(@(x) size(x.target{1}.neural_data.fr,1),single_trial_data);

% get nbr of BDFs (tasks/datasets)
nbr_bdfs                = length(single_trial_data);

% get bin size
bin_size                = single_trial_data{1}.target{1}.bin_size;

% generate new time vector depending on the mode for equalizing trial
% duration that has been chosen  
switch mode
    
    % ---------------------------------------------------------------------
    % resample everything to the duration of the trials of the shortest
    % task
    case 'min_dur'
        new_nbr_bins    = min(nbr_bins_dataset);
        t_new           = 0:bin_size:bin_size*(new_nbr_bins-1);
    
        
    % ---------------------------------------------------------------------
    % resample everything to the duration of the trials of the longest
    % task
    
    case 'max_dur'
        new_nbr_bins    = max(nbr_bins_dataset);
        t_new           = 0 : bin_size : bin_size*(new_nbr_bins-1);

        
    % ---------------------------------------------------------------------
    % cut the trials between time_window(1) and time_window(2)
    
    case 'time_win'
        
        % -----------------------------------------------------------------
        % Some checks to undestand what happens if the code breaks
        
        % check that the time windows are exact multiples of the bin size
        if sum(rem(time_window,bin_size)) > 0
            error('the times in time_window have to be a multiple of the bin size');
        end
        % check that all of the tasks (BDFs) have trials longer than the
        % specified time window
        dur_trials  = bin_size*(nbr_bins_dataset-1);
        if dur_trials < repmat(time_window(2)-time_window(1),1,nbr_bdfs)
            error('>=1 trials are shorter than the specified time window');
        end
        % and that all of the tasks include the time window
        % defined in 'time_window'
        if sum( cellfun( @(x) x.target{1}.t(1) > time_window(1), single_trial_data ) ) > 0
            error('>=1 trials start after time_window(1)');
        end
        if sum( cellfun( @(x) x.target{1}.t(end) < time_window(2), single_trial_data ) ) > 0
            error('>=1 trials finish before time_window(2)');
        end
        
        % -----------------------------------------------------------------
        % define matrix with the indexes of the bins to keep:
        indx_to_keep    = zeros( nbr_bdfs, (time_window(2)-time_window(1)+bin_size)/bin_size );
        start_indx      = cellfun(@(y) find(y,1), cellfun( @(x) x.target{1}.t==time_window(1), ...
                            single_trial_data, 'UniformOutput',false ) );
        end_indx        = cellfun(@(y) find(y,1), cellfun( @(x) x.target{1}.t==time_window(2), ...
                            single_trial_data, 'UniformOutput',false ) );
        
        new_nbr_bins    = numel(start_indx(1):end_indx(1));
                    
        for i = 1:nbr_bdfs
            indx_to_keep(i,:) = start_indx(i):1:end_indx(i);
        end
end



% -------------------------------------------------------------------------
% Retrieve what variables that will be resampled

% get field names, for resampling --- get them from the last target (i.e.
% all the concatenated trials) because some of the 'meta' fields are missing on it 
neural_names            = fieldnames(single_trial_data{1}.target{end}.neural_data);
emg_names               = fieldnames(single_trial_data{1}.target{end}.emg_data);
if isfield(single_trial_data{1}.target{1},'pos')
    pos_names           = fieldnames(single_trial_data{1}.target{end}.pos);
end
if isfield(single_trial_data{1}.target{1},'vel')
    vel_names           = fieldnames(single_trial_data{1}.target{end}.vel);
end
if isfield(single_trial_data{1}.target{1},'force')
    warning('EQUALIZE_SINGLE_TRIAL_DUR: force data length equalization not implemented yet');
end

% get rid of vars that are structs (like 'dim_red' 'neural_data' and in
% 'emg_data')
str_fields_neural       = cellfun(@(x) isstruct(single_trial_data{1}.target{1}.neural_data.(x)), ...
                            neural_names );
str_fields_emg          = cellfun(@(x) isstruct(single_trial_data{1}.target{1}.emg_data.(x)), ...
                            emg_names );

if sum(str_fields_neural) > 0
    % delete the field that is a struct
    neural_names(str_fields_neural) = [];
    neural_dim_red_names = fieldnames(single_trial_data{1}.target{1}.neural_data.dim_red);
end

if sum(str_fields_emg) > 0
    % delete the field that is a struct
    neural_names(str_fields_emg) = [];
    emg_dim_red_names   = fieldnames(single_trial_data{1}.target{1}.emg_data.dim_red);
end



% -------------------------------------------------------------------------
% Resample !


% do for each task
for i = 1:nbr_bdfs
    
    % choose mode for equalizing trial duration
    switch mode
       
        % -----------------------------------------------------------------
        % for the time warping cases
        case {'min_dur','max_dur'}
    
            % create time vector for the original data
            t_orig      = 0 : bin_size : bin_size*(nbr_bins_dataset(i)-1);
    
            
            %--------------------------------------------------------------
            % resample the data for each target
            
            for t = 1:numel(single_trial_data{i}.target)
                
                % -------------------
                % neural data
                for f = 1:numel(neural_names)
                    data_orig   = single_trial_data{i}.target{t}.neural_data.(neural_names{f});
                    % resample if the data are single trials (e.g., not
                    % concatenated), in which case they'll have the same length
                    if size(data_orig,1) == nbr_bins_dataset(i)
                        data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                        single_trial_data{i}.target{t}.neural_data.(neural_names{f}) = ...
                            data_new;
                    end
                end
                
                % -------------------
                % EMG data
                for f = 1:numel(emg_names)
                    data_orig   = single_trial_data{i}.target{t}.emg_data.(emg_names{f});
                    % resample if the data are single trials (e.g., not
                    % concatenated), in which case they'll have the same length
                    if size(data_orig,1) == nbr_bins_dataset(i)
                        data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                        single_trial_data{i}.target{t}.emg_scores.(emg_names{f}) = ...
                            data_new;
                    end
                end
                
                % -------------------
                % dim reduced neural data
                if exist('neural_dim_red_names','var')
                    for f = 1:numel(neural_dim_red_names)
                        data_orig = single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                            single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f}) = ...
                                data_new;
                        end
                    end
                end
                
                % -------------------
                % dim reduced EMG data
                if exist('emg_dim_red_names','var')
                    for f = 1:numel(emg_dim_red_names)
                        data_orig = single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                            single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{f}) = ...
                                data_new;
                        end
                    end
                end
                
                % -------------------        
                % Pos data
                if exist('pos_names','var')
                    for f = 1:numel(pos_names)
                        data_orig   = single_trial_data{i}.target{t}.pos.(pos_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                            single_trial_data{i}.target{t}.pos.(pos_names{f}) = ...
                                data_new;
                        end
                    end
                end
                
                % -------------------
                % Vel data
                if exist('vel_names','var')
                    for f = 1:numel(vel_names)
                        data_orig   = single_trial_data{i}.target{t}.vel.(vel_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                            single_trial_data{i}.target{t}.pos.(vel_names{f}) = ...
                                data_new;
                        end
                    end
                end
            end
            
            
        % -----------------------------------------------------------------
        % for cutting the trials between certain t_init and t_end
        
        case 'time_win'
            
            %--------------------------------------------------------------
            % resample the data for each target
            for t = 1:numel(single_trial_data{t})
               
                % cut the time vector
                single_trial_data{i}.target{t}.t = single_trial_data{i}.target{t}.t(indx_to_keep(i,:));
                
                % -------------------
                % neural data
                for f = 1:numel(neural_names)
                    data_orig       = single_trial_data{i}.target{t}.neural_data.(neural_names{f});
                    % resample if the data are single trials (e.g., not
                    % concatenated), in which case they'll have the same length
                    if size(data_orig,1) == nbr_bins_dataset(i)
                        % choose the selected bins (note that matlab is
                        % smart and if the matrix is 2D this code works,
                        % but also works with the 3D FR matrices
                        data_new    = single_trial_data{i}.target{t}.neural_data.(neural_names{f})(indx_to_keep(i,:),:,:);
                        single_trial_data{i}.target{t}.neural_data.(neural_names{f}) = data_new;
                    end     
                end
                
                % -------------------
                % EMG data
                for f = 1:numel(emg_names)
                    data_orig       = single_trial_data{i}.target{t}.emg_data.(emg_names{f});
                    % resample if the data are single trials (e.g., not
                    % concatenated), in which case they'll have the same length
                    if size(data_orig,1) == nbr_bins_dataset(i)
                        data_new    = single_trial_data{i}.target{t}.emg_data.(emg_names{f})(indx_to_keep(i,:),:,:);
                        single_trial_data{i}.target{t}.emg_data.(emg_names{f}) = data_new;
                    end
                end
                
                % -------------------
                % dim reduced neural data
                if exist('neural_dim_red_names','var')
                    for f = 1:numel(neural_dim_red_names)
                        data_orig   = single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new = single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f})...
                                            (indx_to_keep(i,:),:,:);
                            single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f}) = ...
                                data_new;
                        end
                    end
                end
                
                % -------------------
                % dim reduced EMG data
                if exist('emg_dim_red_names','var')
                    for f = 1:numel(emg_dim_red_names)
                        data_orig   = single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new = single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{f})...
                                            (indx_to_keep(i,:),:,:);
                            single_trial_data{i}.target{t}.emg_data.dim_red.(emg_dim_red_names{f}) = ...
                                data_new;
                        end
                    end
                end
                
                % -------------------        
                % Pos data
                if exist('pos_names','var')
                    for f = 1:numel(pos_names)
                        data_orig   = single_trial_data{i}.target{t}.pos.(pos_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new = single_trial_data{i}.target{t}.pos.(pos_names{f})(indx_to_keep(i,:),:);
                            single_trial_data{i}.target{t}.pos.(pos_names{f}) = data_new;
                        end
                    end
                end
                
                % -------------------        
                % Vel data
                if exist('vel_names','var')
                    for f = 1:numel(vel_names)
                        data_orig   = single_trial_data{i}.target{t}.vel.(vel_names{f});
                        % resample if the data are single trials (e.g., not
                        % concatenated), in which case they'll have the same length
                        if size(data_orig,1) == nbr_bins_dataset(i)
                            data_new = single_trial_data{i}.target{t}.vel.(vel_names{f})(indx_to_keep(i,:),:);
                            single_trial_data{i}.target{t}.vel.(vel_names{f}) = data_new;
                        end
                    end
                end
            end
    end
end

