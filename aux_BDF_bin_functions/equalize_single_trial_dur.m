% 
% Take an array of target averaged responses for several tasks obtained
% with single_trial_analysis_dim_red and make them all have the same
% duration (number of bins)
%
%   function single_trial_data   = equalize_single_trial_dur( single_trial_data )
%
%

function single_trial_data   = equalize_single_trial_dur( single_trial_data, varargin )


% -------------------------------------------------------------------------
% Decide to what number of bins to resample all trial averaged data


% get number of samples of the average responses for each dataset
bins_dataset        = cellfun(@(x) size(x.target{1}.neural_data.fr,1),single_trial_data);

% resample everything to the duration of the longest
new_nbr_bins        = max(bins_dataset);
t_new               = 0:single_trial_data{1}.target{1}.bin_size:...
                        single_trial_data{1}.target{1}.bin_size*(new_nbr_bins-1);
            
   
% -------------------------------------------------------------------------
% Retrieve the variables that will be resampled

% get field names, for resampling --- get them from the last target (i.e.
% all the concatenated trials) because some of the 'meta' fields are missing on it 
neural_names        = fieldnames(single_trial_data{1}.target{end}.neural_data);
emg_names           = fieldnames(single_trial_data{1}.target{end}.emg_data);
if isfield(single_trial_data{1}.target{1},'pos')
    pos_names       = fieldnames(single_trial_data{1}.target{end}.pos);
end
if isfield(single_trial_data{1}.target{1},'vel')
    vel_names       = fieldnames(single_trial_data{1}.target{end}.vel);
end
if isfield(single_trial_data{1}.target{1},'force')
    warning('EQUALIZE_SINGLE_TRIAL_DUR: force data length equalization not implemented yet');
end

% get rid of vars that are structs (like 'dim_red' 'neural_data' and in
% 'emg_data')
str_fields_neural   = cellfun(@(x) isstruct(single_trial_data{1}.target{1}.neural_data.(x)), neural_names );
str_fields_emg      = cellfun(@(x) isstruct(single_trial_data{1}.target{1}.emg_data.(x)), emg_names );

if sum(str_fields_neural) > 0
    % delete the field that is a struct
    neural_names(str_fields_neural) = [];
    neural_dim_red_names = fieldnames(single_trial_data{1}.target{1}.neural_data.dim_red);
end

if sum(str_fields_emg) > 0
    % delete the field that is a struct
    neural_names(str_fields_emg) = [];
    emg_dim_red_names = fieldnames(single_trial_data{1}.target{1}.emg_data.dim_red);
end


% -------------------------------------------------------------------------
% Resample !


% do for each task
for i = 1:numel(single_trial_data)
    
    % create time vector for the original data
    t_orig          = 0:single_trial_data{i}.target{1}.bin_size:...
                        single_trial_data{i}.target{1}.bin_size*(bins_dataset(i)-1);
    
    %----------------------------------
 	% resample the data for each target
    for t = 1:numel(single_trial_data{i}.target)
        % -------------------
        % for the neural data
        for f = 1:numel(neural_names)
            data_orig   = single_trial_data{i}.target{t}.neural_data.(neural_names{f});
            % resample if the data are single trials (e.g., not
            % concatenated), in which case they'll have the same length
            if size(data_orig,1) == length(t_orig)
                data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                single_trial_data{i}.target{t}.neural_data.(neural_names{f}) = ...
                    data_new;
            end
        end
        % -------------------
        % for the EMG data
        for f = 1:numel(emg_names)
            data_orig   = single_trial_data{i}.target{t}.emg_data.(emg_names{f});
            % resample if the data are single trials (e.g., not
            % concatenated), in which case they'll have the same length
            if size(data_orig,1) == length(t_orig)
                data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                single_trial_data{i}.target{t}.emg_scores.(emg_names{f}) = ...
                    data_new;
            end
        end
        % -------------------
        % for the dim reduced neural data
        if exist('neural_dim_red_names','var')
            for f = 1:numel(neural_dim_red_names)
                data_orig = single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f});
                % resample if the data are single trials (e.g., not
                % concatenated), in which case they'll have the same length
                if size(data_orig,1) == length(t_orig)
                    data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                    single_trial_data{i}.target{t}.neural_data.dim_red.(neural_dim_red_names{f}) = ...
                        data_new;
                end
            end
        end
        % -------------------
        % for the dim reduced EMG data
        if exist('emg_dim_red_names','var')
            for f = 1:numel(emg_dim_red_names)
                data_orig = single_trial_data{i}.target{t}.neural_data.dim_red.(emg_dim_red_names{f});
                % resample if the data are single trials (e.g., not
                % concatenated), in which case they'll have the same length
                if size(data_orig,1) == length(t_orig)
                    data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                    single_trial_data{i}.target{t}.neural_data.dim_red.(emg_dim_red_names{f}) = ...
                        data_new;
                end
            end
        end
        % -------------------        
        % for the Pos data
        if exist('pos_names','var')
            for f = 1:numel(pos_names)
                data_orig   = single_trial_data{i}.target{t}.pos.(pos_names{f});
                % resample if the data are single trials (e.g., not
                % concatenated), in which case they'll have the same length
                if size(data_orig,1) == length(t_orig)
                    data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                    single_trial_data{i}.target{t}.pos.(pos_names{f}) = ...
                        data_new;
                end
            end
        end
        % -------------------
        % for the Vel data
        if exist('vel_names','var')
            for f = 1:numel(vel_names)
                data_orig   = single_trial_data{i}.target{t}.vel.(vel_names{f});
                % resample if the data are single trials (e.g., not
                % concatenated), in which case they'll have the same length
                if size(data_orig,1) == length(t_orig)
                    data_new    = interp1(t_orig,data_orig,t_new,'linear','extrap');
                    single_trial_data{i}.target{t}.pos.(vel_names{f}) = ...
                        data_new;
                end
            end
        end
    end
end

