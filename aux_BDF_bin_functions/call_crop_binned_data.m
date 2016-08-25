%
% Crop binned_data from a 'Wrist Flexion', 'Multigadget' or 'Ball' task in
% intervals defined by two words. The function can take a BDF that will be
% first converted into a binned_data struct
%
%   function cropped_binned_data = call_crop_binned_data( data_struct, ...
%                                   word_i, word_f, task, varargin )
%
% Inputs (opt)              : [default]
%   data_struct             : array of BDFs or binned_data structs
%   word_i                  : word that defines the beginning of the
%                               interval ('start','ot_on','go')
%   word_f                  : word that defines the end of the interval
%                               ('ot_on','go','end','R')
%   task                    : the behavior: 'wf': wrist flexion, 'mg':
%                               multi-gadget, 'bd': ball task
%   (bin_size)              : [0.05 s] bin size for binning a BDF
%
% Outputs
%   cropped_binned_data     : cropped binned_data struct
%
%
% Usage:
%   cropped_binned_data = call_crop_binned_data( binned_data, word_i, word_f, task )
%   cropped_binned_data = call_crop_binned_data( bdf, word_i, word_f, task )
%   cropped_binned_data = call_crop_binned_data( bdf, word_i, word_f, task, bin_size )
  

function cropped_binned_data = call_crop_binned_data( data_struct, word_i, word_f, task, varargin )


% see if the data_struct is of type bdf or binned_data. If it is a BDF,
% convert it to a binned_data struct
if ~isfield(data_struct,'timeframe')
    % get desired bin size
    if nargin == 5
        bin_pars.binsize    = varargin{1};
    % or set it to default
    else
        bin_pars.binsize    = 0.05;
    end
    
    % bin each BDF
    for i = 1:length(data_struct)
        % start and stop times are set to ensure compatibility with the
        % dim_reduction code
        bin_pars.starttime  = ceil(data_struct(i).pos(1,1)/bin_pars.binsize)*...
                                    bin_pars.binsize;
        bin_pars.stoptime   = floor(data_struct(i).pos(end,1)/bin_pars.binsize)*...
                                    bin_pars.binsize;
        binned_data_array(i) = convertBDF2binned(data_struct(i),bin_pars);
    end
else
    binned_data_array       = data_struct;
end

clear data_struct;
    

% get nbr of binned_data structs
nbr_bdfs                    = length(binned_data_array);


for i = 1:nbr_bdfs

    % get trial table
    trial_table             = binned_data_array.trialtable;
    
    % get the column of the words. Note that the organization of the trial
    % table depends on the task
    switch task{i}
        % WRIST FLEXION
        case {'wf','iso8','iso','wm','spr'}
            % first word for cropping
            switch word_i
                case 'start'
                    indx_i  = 1;
                case 'ot_on'
                    indx_i  = 6;
                case 'go'
                    indx_i  = 7;
                otherwise
                    error([word_i ' not supported for ' task ' task']);
            end
            % last word for cropping
            switch word_f
                case 'ot_on'
                    indx_f  = 6;
                case 'go'
                    indx_f  = 7;
                case 'end'
                    indx_f  = 8;
                case 'R'
                    indx_f  = 8; % code will then look at whether the monkey got a reward
                otherwise
                    error([word_f ' not supported for this task']);
            end
        % MULTIGADGET 
        case {'mg','mg-pt'}
            
%             % get nbr cols of the trial table --in very old datasets there
%             % were only 7 columns, where 1 seems to correspond to 'start', 2
%             % to 'ot_on', 3 to 'touchpad off' (?), 4 to target number, 5 to
%             % contact sensor time (picking up ball) (?), 6 to trial end,
%             % and 7 the trial end word
%             cols_trial_table_mg = size(trial_table,2);
            
            % first word for cropping
            switch word_i
                case 'start'
                    indx_i  = 1;
                case 'ot_on'
                    indx_i  = 2;
                    warning('using touchpad time as OT ON');
                case 'go'
                    indx_i  = 3;
                otherwise
                    error([word_i ' not supported for ' task ' task']);
            end
            % last word for cropping
            switch word_f
                case 'ot_on'
                    indx_f  = 2;
                case 'go'
                    indx_f  = 3;
                case 'end'
%                     % for compatibility with really old files
%                     if cols_trial_table_mg == 12
%                         indx_f  = 11;
%                     elseif cols_trial_table_mg == 7
%                         indx_f  = 6;
%                     end
                    indx_f  = 11;
                case 'R'
%                     % for compatibility with really old files                    
%                     if cols_trial_table_mg == 12
%                         indx_f  = 11; % code will then look at whether the monkey got a reward
%                     elseif cols_trial_table_mg == 7
%                         indx_f  = 6;
%                     end
                    indx_f  = 11;
                otherwise
                    error([word_f ' not supported for this task']);
            end
        % BALL DEVICE
        case 'ball'
            % first word for cropping
            switch word_i
                case 'start'
                    indx_i  = 1;
                case 'ot_on'
                    indx_i  = 2;
                    warning('using touchpad time as OT ON');
                case 'go'
                    indx_i  = 3;
                otherwise
                    error([word_i ' not supported for ' task ' task']);
            end
            % last word for cropping
            switch word_f
                case 'ot_on'
                    indx_f  = 2;
                case 'go'
                    indx_f  = 3;
                case 'end'
                    indx_f  = 6;
                case 'R'
                    indx_f  = 6; % code will then look at whether the monkey got a reward
                otherwise
                    error([word_f ' not supported for this task']);
            end
    end

    % get trial table and store it in a N x 2 matrix with times for
    % cropping 
    cropping_times          = [trial_table(:,indx_i), trial_table(:,indx_f)];

    % if the end word is 'R' (reward), get rid of the trials without a reward
    if word_f == 'R'
        switch task{i}
            case {'wf','iso8','iso','wm','spr'}
                col_R       = 9;
            case {'mg','mg-pt'}
                % for compatibility with really old files
%                 switch cols_trial_table_mg
%                     case 12
%                         col_R = 12;
%                     case 7
%                         col_R = 7;
%                 end
                col_R       = 12;
            case 'ball'
                col_R       = 7;
        end     
        cropping_times(trial_table(:,col_R) ~= double('R'),:) = [];
        binned_data_array(i).trialtable( binned_data_array(i).trialtable(:,col_R) ...
            ~= double('R'), : ) = [];
    end
    
    % there's an error in our behavior, which sometimes makes the
    % trial_start time == -1
    if ~isempty(find(cropping_times(:,1)==-1,1))
       cropping_times(cropping_times(:,1)==-1,:) = [];
    end

    % call cropping function
    cropped_binned_data(i)  = crop_binned_data( binned_data_array(i), cropping_times );

    
    clear trial_table cropping_times;
end