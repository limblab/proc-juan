%
% Find cutting times for binned data based on the trial table.
%
%    [cut_b, cut_t] = find_cutting_times( trial_table, varargin )
%
%
% Inputs (opt)      : [default]
%   trial_table     : the task's trial table (e.g., from a binned_data
%                       file)
%   ('task')        : ['wf'] the task ('wf','iso','spr','mg','mg-pt','ball')
%   ('word_i')      : ['ot_on'] data window start word ('start','ot_on','go') 
%   ('word_f')      : ['R'] data window end word ('ot_on','go','end','R')
%   ('bin_size')    : [0.02] bin size (ms)
%
%


function [cut_b, cut_t] = find_cutting_times( trial_table, varargin )


% default parameters
params                  = struct('task',        'wf',...
                                'word_i',       'ot_on',...
                                'word_f',       'R',...
                                'bin_size',     0.02);
                    

% read input parameters and replace default values where necessary
param_names             = fieldnames(params);
for p = reshape(varargin,2,[])
    if any(strcmp(p{1},param_names))
        params.(p{1})   = p{2};
    else
        error([p{1} ' is not a recognized parameter name']);
    end
end



% 1. retrieve columns in the trial table that have the times of the desired
% words

switch params.task
    
    %--------------------------
    % WRIST FLEXION
    case {'wf','iso8','iso','wm','spr'}
        % first word for cropping
        switch params.word_i
            case 'start'
                indx_i  = 1;
            case 'ot_on'
                indx_i  = 6;
            case 'go'
                indx_i  = 7;
            otherwise
                error([params.word_i ' not supported for ' task ' task']);
        end
        % last word for cropping
        switch params.word_f
            case 'ot_on'
                indx_f  = 6;
            case 'go'
                indx_f  = 7;
            case 'end'
                indx_f  = 8;
            case 'R'
                indx_f  = 8; % code will then look at whether the monkey got a reward
            otherwise
                error([params.word_f ' not supported for this task']);
        end
        
    %--------------------------
    % MULTIGADGET
    case {'mg','mg-pt'}
        
        %             % get nbr cols of the trial table --in very old datasets there
        %             % were only 7 columns, where 1 seems to correspond to 'start', 2
        %             % to 'ot_on', 3 to 'touchpad off' (?), 4 to target number, 5 to
        %             % contact sensor time (picking up ball) (?), 6 to trial end,
        %             % and 7 the trial end word
        %             cols_trial_table_mg = size(trial_table,2);
        
        % first word for cropping
        switch params.word_i
            case 'start'
                indx_i  = 1;
            case 'ot_on'
                indx_i  = 2;
                warning('using touchpad time as OT ON');
            case 'go'
                indx_i  = 3;
            otherwise
                error([params.word_i ' not supported for ' task ' task']);
        end
        % last word for cropping
        switch params.word_f
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
                error([params.word_f ' not supported for this task']);
        end
        
    %--------------------------
    % BALL DEVICE
    case 'ball'
        % first word for cropping
        switch params.word_i
            case 'start'
                indx_i  = 1;
            case 'ot_on'
                indx_i  = 2;
                warning('using touchpad time as OT ON');
            case 'go'
                indx_i  = 3;
            otherwise
                error([params.word_i ' not supported for ' task ' task']);
        end
        % last word for cropping
        switch params.word_f
            case 'ot_on'
                indx_f  = 2;
            case 'go'
                indx_f  = 3;
            case 'end'
                indx_f  = 6;
            case 'R'
                indx_f  = 6; % code will then look at whether the monkey got a reward
            otherwise
                error([params.word_f ' not supported for this task']);
        end
end


% 2. create a (trials x 2) matrix that has the indexes for cutting
cut_t           = [trial_table(:,indx_i), trial_table(:,indx_f)];


% 2.b. if the end word is 'R' (reward), get rid of the trials without a
% reward 
if params.word_f == 'R'
    switch params.task
        case {'wf','iso8','iso','wm','spr'}
            col_R       = 9;
        case {'mg','mg-pt'}
%            % for compatibility with really old files
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
    
    cut_t(trial_table(:,col_R) ~= double('R'),:) = [];
    
    % there's an error in our behavior, which sometimes makes the
    % trial_start time == -1
    if ~isempty(find(cut_t(:,1)==-1,1))
       cut_t(cut_t(:,1)==-1,:) = [];
    end
    
    % turn cut times into bin number times
    cut_b               = zeros(size(cut_t));
    cut_b(:,1)          = ceil(cut_t(:,1)/params.bin_size);
    cut_b(:,2)          = floor(cut_t(:,2)/params.bin_size);
end
