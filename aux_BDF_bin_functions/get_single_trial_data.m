%
% get single trial data from binned_data struct.
%
%   single_trial_data = get_single_trial_data( binned_data, task, varargin )
%
%
% Inputs (opt)      : [default]
%   binned_data     : binned_data struct
%   neural_chs      : neural chs that will be used in the analysis. All if
%                       empty
%   task            : task the monkey performed. So far compatible with
%                       'wf', 'mg', 'ball'
%   (trial_norm)    : ['stretch'] how the trials will be averaged:
%                       'stretch': resampled ot the duration of the
%                       shortest trial; 'min_dur': cropped to the duration
%                       of the shortest trial
%   (neural_chs)    : [all] Neural chs that will be used. All if empty
%   (emg_chs)       : [all] EMG chs that will be used. All if empty
%   (w_i)           : ['ot_on'] word that defines trial start ('ot_on', 
%                       'go', 'start')
%   (w_f)           : ['R'] word that defines trial end ('ot_on', 'go', 
%                       'end', 'R')
%   (dim_red_FR)    : [~] dim_red_FR struct
%   (dim_red_emg)   : [~] dim_red_emg struct
%
%
% Outputs:
%   single_trial_data : struct with single trial data, organized by target
%                       (one target per field of the cell array, plus a
%                       last field with all the targets)
%
%
%
% function single_trial_data = get_single_trial_data( binned_data, task )
% function single_trial_data = get_single_trial_data( binned_data, task, trial_norm )
% function single_trial_data = get_single_trial_data( binned_data, task, neural_chs, emg_chs )
% function single_trial_data = get_single_trial_data( binned_data, task, neural_chs, emg_chs, trial_norm )
% function single_trial_data = get_single_trial_data( binned_data, task, neural_chs, emg_chs, w_i, w_f )
% function single_trial_data = get_single_trial_data( binned_data, task, neural_chs, emg_chs, trial_norm, w_i, w_f )
%
% function single_trial_data = get_single_trial_data( binned_data, task, neural_chs, emg_chs, trial_norm, w_i, w_f, dim_red_FR, dim_red_emg )
%
%
% ToDo's
%   - allow the user to pass a BDF instead of a binned_data struct, since
%   the functions inside can take binned_data structs
%


function single_trial_data = get_single_trial_data( binned_data, task, varargin )


% -------------------------------------------------------------------------
% read input parameters

if nargin == 3
    trial_norm              = varargin{1};
elseif nargin == 4
    neural_chs              = varargin{1};
    emg_chs                 = varargin{2};
elseif nargin == 5
    neural_chs              = varargin{1};
    emg_chs                 = varargin{2};
    trial_norm              = varargin{3};
elseif nargin == 6    
    neural_chs              = varargin{1};
    emg_chs                 = varargin{2};
    w_i                     = varargin{3};
    w_f                     = varargin{4};
end
if nargin >= 7
    neural_chs              = varargin{1};
    emg_chs                 = varargin{2};
    trial_norm              = varargin{3};
    w_i                     = varargin{4};
    w_f                     = varargin{5};
end
if nargin == 9
    dim_red_FR              = varargin{6};
    dim_red_emg             = varargin{7};
    dim_red_data_yn         = true; % flag for processing these data
end

% set defaults for not passed params
if ~exist('w_i','var')
    w_i                     = 'ot_on';
    w_f                     = 'R';
end
if ~exist('trial_norm','var')
    trial_norm              = 'stretch';
end
if ~exist('neural_chs','var')
   neural_chs               = 1:size(binned_data.spikeratedata,2); 
end
if ~exist('emg_chs','var')
    emg_chs                 = 1:size(binned_data.emgdatabin,2);
end
if ~exist('dim_red_FR','var')
    dim_red_data_yn = false;
end

% fill neural_chs and emg_chs to all chs if empty
if isempty(neural_chs)
   neural_chs               = 1:size(binned_data.spikeratedata,2); 
end
if isempty(emg_chs)
    emg_chs                 = 1:size(binned_data.emgdatabin,2);
end


% -------------------------------------------------------------------------
% only keep the neural/emg channels we want


nbr_neural_chs              = length(neural_chs);
nbr_emgs                    = length(emg_chs);

binned_data.spikeratedata   = binned_data.spikeratedata(:,neural_chs);
if isfield(binned_data,'smoothedspikerate')
    binned_data.smoothedspikerate = binned_data.smoothedspikerate(:,neural_chs);
end
binned_data.emgdatabin      = binned_data.emgdatabin(:,emg_chs);


% -------------------------------------------------------------------------
% chop continuous data to task-related data only


cropped_binned_data         = call_crop_binned_data( binned_data, w_i, w_f, task );

% get bin size
bin_size                    = min(diff(binned_data.timeframe));
bin_size                    = round(bin_size*100)/100;

% get information about targets: how many and coordinates
[nbr_targets, target_coord, target_id] = get_targets( cropped_binned_data, task, false );

% chop the dim_red_data to the same bins, if present
if dim_red_data_yn
    [dim_red_FR, dim_red_emg] = crop_dim_red_data_based_on_binned( cropped_binned_data, ...
                                dim_red_FR, dim_red_emg );
end


% ------------------------------------------------------------------------
% Get trial info


% retrieve the intervals that correspond to each trial
diff_timeframe              = diff(cropped_binned_data.timeframe);

% get the end of each trial
end_indx                    = find(diff_timeframe >= 2*bin_size);
% add the end of tehe last trial
end_indx                    = [end_indx; length(cropped_binned_data.timeframe)];

% get the beginning of each trial
start_indx                  = [1; end_indx(1:end-1)+1];
% get the number of trials
nbr_trials                  = length(end_indx);
% get trial duration (for all trials)
trial_dur                   = zeros(nbr_trials,1);
trial_dur(1)                = end_indx(1)*bin_size;
trial_dur(2:nbr_trials)     = diff(end_indx)*bin_size;


% ToDo: improve this
% check that the trials have been retrieved as they should
% Note: because of the binning, sometimes the last target(s) are missing
while cropped_binned_data.trialtable(end,end-1)-cropped_binned_data.timeframe(end) > 0
    warning('last target fall oustide the binned data')
    target_id(end)          = [];
    cropped_binned_data.trialtable(end,:) = [];
end


% ------------------------------------------------------------------------
% Create struct with trial-related data

% init cell array with a field per target + an additional field with all of
% the targets concatenated
STD                         = cell(1,nbr_targets+1);


% get duration to which trials will be cut
switch trial_norm
    % only take the data from t = 0 : min(duration_all_trials_this_target)
    case 'min_dur'
        st_dur              = min(trial_dur);
    % 'stretch' all trials in the time domain, so they last the same as
    % the shortest trial
    case 'stretch'
        % st_dur      = ceil(max((trial_dur(trials_this_target)))/bin_size)*bin_size;
        st_dur              = ceil(min(trial_dur)/bin_size)*bin_size;
end


% FIX: in some old files there is one target that shouldn't have been
% there, and that only appears one time
% --I've implemented a fix, in which the id of that target in target_id =
% NaN, so it can be corrected when generating the variable that points to
% the bins in cropped_binned_data that correspond to each trial
wrong_target_indx           = find(isnan(target_id));
if length(wrong_target_indx) > 1
    error('get_single_trial_data.m: more than one wrong target not supported yet');
end



% ---------------------------------
% get data for each target

for i = 1:nbr_targets
    
    % get the trials that correspond to this target (target id = 0,1,2...)
    trials_this_tgt         = find(target_id == i-1);
    nbr_trials_this_tgt     = length(trials_this_tgt);
    
    % get nbr of bins per trial
    nbr_bins                = round(st_dur/bin_size);

    bin_indx_p_trial        = zeros(nbr_bins,length(trials_this_tgt));
    
    
    % FIX: create an offset variable that will be adapted if trials are
    % excluded from the analysis, as it gives problems when creating
    % bin_indx. This is necessary because in some old data there are
    % targets that appear just once and that had not been defined 
    trial_offset            = 0;

    
    % ---------------------------------
    % preallocate data matrices
        
    % for the raw data
    STD{i}.neural_data.fr   = zeros( nbr_bins,...
                                nbr_neural_chs, nbr_trials_this_tgt );
    if isfield(cropped_binned_data,'smoothedspikerate')
        STD{i}.neural_data.smoothed_fr = zeros( nbr_bins,...
                                nbr_neural_chs, nbr_trials_this_tgt );
    end
    
    STD{i}.emg_data.emg     = zeros( nbr_bins,...
                                nbr_emgs, nbr_trials_this_tgt );
	if ~isempty(cropped_binned_data.cursorposbin)
        STD{i}.pos.data     = zeros( nbr_bins,...
                                2, nbr_trials_this_tgt );
        STD{i}.vel.data     = zeros( nbr_bins,...
                                3, nbr_trials_this_tgt );
    end
    
    % some statistics (mean, SD)
    STD{i}.neural_data.mn 	= zeros( nbr_bins,...
                                nbr_neural_chs );
    if isfield(cropped_binned_data,'smoothedspikerate')
        STD{i}.neural_data.smoothed_fr_mn = zeros( nbr_bins,...
                                nbr_neural_chs );
    end
    STD{i}.emg_data.mn      = zeros( nbr_bins,...
                                nbr_emgs );
    if ~isempty(cropped_binned_data.cursorposbin)
        STD{i}.pos.mn       = zeros( nbr_bins, 2 );
        STD{i}.vel.mn       = zeros( nbr_bins, 3 );
    end
    
    STD{i}.neural_data.sd   = zeros( nbr_bins,...
                                nbr_neural_chs );
	if isfield(cropped_binned_data,'smoothedspikerate')
        STD{i}.neural_data.smoothed_fr_sd = zeros( nbr_bins,...
                                nbr_neural_chs );
    end
    STD{i}.emg_data.sd      = zeros( nbr_bins,...
                                nbr_emgs );
    if ~isempty(cropped_binned_data.cursorposbin)            
        STD{i}.pos.sd       = zeros( nbr_bins, 2 );
        STD{i}.vel.sd       = zeros( nbr_bins, 3 );
    end
    
    % for the dimensionality-reduced data
    if dim_red_data_yn
       
        STD{i}.neural_data.dim_red.scores = zeros( nbr_bins,...
                                nbr_neural_chs, nbr_trials_this_tgt );
        STD{i}.neural_data.dim_red.mn = zeros( nbr_bins,...
                                nbr_neural_chs );
        STD{i}.neural_data.dim_red.sd = zeros( nbr_bins,...
                                nbr_neural_chs );
                                
        STD{i}.emg_data.dim_red.scores = zeros( nbr_bins,...
                                nbr_emgs, nbr_trials_this_tgt );
        STD{i}.emg_data.dim_red.mn = zeros( nbr_bins,...
                                nbr_emgs );
        STD{i}.emg_data.dim_red.sd = zeros( nbr_bins,...
                                nbr_emgs );
    end
    
    % for a matrix that will include the bin number after cropping
    STD{i}.bin_indx         = zeros(1,nbr_trials_this_tgt*nbr_bins);
    
    
    % ---------------------------------
    % fill with values
    
    for t = 1:nbr_trials_this_tgt
        % find trial start
        aux_start           = start_indx(trials_this_tgt(t));
        
        switch trial_norm 
            
            case 'min_dur'
                % find trial end: for 'min_dur' happens at the end of the
                % shortest trial for this target
                aux_end     = aux_start + nbr_bins - 1;
                
                % raw data
                STD{i}.neural_data.fr(:,:,t) = cropped_binned_data.spikeratedata(...
                                aux_start:aux_end,:);
                if isfield(cropped_binned_data,'smoothedspikerate')
                    STD{i}.neural_data.smoothed_fr(:,:,t) = cropped_binned_data.smoothedspikerate(...
                                aux_start:aux_end,:);                    
                end
                
                STD{i}.emg_data.emg(:,:,t) = cropped_binned_data.emgdatabin(...
                                aux_start:aux_end,:);
                
                if ~isempty(cropped_binned_data.cursorposbin)
                    STD{i}.pos.data(:,:,t) = ...
                        cropped_binned_data.cursorposbin(aux_start:aux_end,:);
                    STD{i}.vel.data(:,:,t) = ...
                        cropped_binned_data.velocbin(aux_start:aux_end,:);
                end
                
                % dimensionality reduced data
                if dim_red_data_yn
                    STD{i}.neural_data.dim_red.scores(:,:,t) = dim_red_FR.scores(...
                                aux_start:aux_end,:);
                    STD{i}.emg_data.dim_red.scores(:,:,t) = dim_red_emg.scores(...
                                aux_start:aux_end,:);
                end
                                
            case 'stretch'
                
                % -----------------
                % 1. find trial end: for 'stretch is the real end of the
                % trial; afterwards the code will make all the trials equal
                % length
                aux_end             = aux_start + trial_dur(trials_this_tgt(t))/bin_size - 1;
                
                % -----------------
                % 2. for interpolating, define time axis of the data
                t_orig              = 1:trial_dur(trials_this_tgt(t))/bin_size;
                t_new               = 1:nbr_bins;
                
                fr_orig             = cropped_binned_data.spikeratedata(aux_start:aux_end,:);
                if isfield(cropped_binned_data,'smoothedspikerate')
                    smoothed_fr_orig = cropped_binned_data.smoothedspikerate(aux_start:aux_end,:);
                end
                emg_orig            = cropped_binned_data.emgdatabin(aux_start:aux_end,:);
                if ~isempty(cropped_binned_data.cursorposbin)
                    pos_data_orig   = cropped_binned_data.cursorposbin(aux_start:aux_end,:);
                    vel_data_orig   = cropped_binned_data.velocbin(aux_start:aux_end,:);
                end
                
                if dim_red_data_yn
                    neural_scores_orig = dim_red_FR.scores(...
                                aux_start:aux_end,:);
                    emg_scores_orig = dim_red_emg.scores(...
                                aux_start:aux_end,:);
                end
                
                % -----------------
                % 3. and now interpolate
                fr_new              = interp1(t_orig,fr_orig,t_new,'linear','extrap');
                if isfield(cropped_binned_data,'smoothedspikerate')
                    smoothed_fr_new = interp1(t_orig,smoothed_fr_orig,t_new,'linear','extrap');
                end
                emg_new             = interp1(t_orig,emg_orig,t_new,'linear','extrap');
                if ~isempty(cropped_binned_data.cursorposbin)
                    pos_data_new    = interp1(t_orig,pos_data_orig,t_new,'linear','extrap');
                    vel_data_new    = interp1(t_orig,vel_data_orig,t_new,'linear','extrap');
                end

                if dim_red_data_yn
                    neural_scores_new = interp1(t_orig,neural_scores_orig,t_new,'linear','extrap');
                    emg_scores_new  = interp1(t_orig,emg_scores_orig,t_new,'linear','extrap');
                end                
                
                % -----------------                
                % 4. store values in single_trial_data struct
                STD{i}.neural_data.fr(:,:,t) = fr_new;
                if isfield(cropped_binned_data,'smoothedspikerate')
                    STD{i}.neural_data.smoothed_fr(:,:,t) = smoothed_fr_new;
                end
                STD{i}.emg_data.emg(:,:,t) = emg_new;
                if ~isempty(cropped_binned_data.cursorposbin)
                    STD{i}.pos.data(:,:,t) = pos_data_new;
                    STD{i}.vel.data(:,:,t) = vel_data_new;
                end
                
                if dim_red_data_yn
                    STD{i}.neural_data.dim_red.scores(:,:,t) = neural_scores_new;
                    STD{i}.emg_data.dim_red.scores(:,:,t) = emg_scores_new;
                end
        end
        
        % FIX: check if we need to add an offset to the trial counter
        if ~isempty(wrong_target_indx)
            if trials_this_tgt(t) == trials_this_tgt( find(trials_this_tgt > wrong_target_indx,1) )
                trial_offset        = trial_offset + 1;
            end
        end
        
        % bins (of the cropped_binned_data struct) that correspond
        % to this trial
        % -- FIX: I had to subtract the trial offset to deal with the
        % problems generated by undesired targets
        STD{i}.bin_indx( (t-1)*nbr_bins+1 : t*nbr_bins ) = ...
                    (trials_this_tgt(t)-1-trial_offset)*nbr_bins+1 : ...
                    (trials_this_tgt(t)-trial_offset)*nbr_bins;
                
        % add the bins that correspond to this trial to bins_p_trial
        STD{i}.bin_indx_p_trial(:,t) = (trials_this_tgt(t)-1-trial_offset)*nbr_bins+1 : ...
                    (trials_this_tgt(t)-trial_offset)*nbr_bins;
    end
    
    % -----------------    
    % calculate mean and SD
    
    STD{i}.neural_data.mn   = mean(STD{i}.neural_data.fr,3);
    STD{i}.neural_data.sd   = std(STD{i}.neural_data.fr,0,3);
    
    if isfield(cropped_binned_data,'smoothedspikerate')
        STD{i}.neural_data.smoothed_fr_mn = mean(STD{i}.neural_data.smoothed_fr,3);
        STD{i}.neural_data.smoothed_fr_sd = std(STD{i}.neural_data.smoothed_fr,0,3);
    end
    
    STD{i}.emg_data.mn      = mean(STD{i}.emg_data.emg,3);
    STD{i}.emg_data.sd      = std(STD{i}.emg_data.emg,0,3);
    
    if ~isempty(cropped_binned_data.cursorposbin)
        STD{i}.pos.mn      	= mean(STD{i}.pos.data,3);
        STD{i}.pos.sd       = std(STD{i}.pos.data,0,3);

        STD{i}.vel.mn       = mean(STD{i}.vel.data,3);
        STD{i}.vel.sd       = std(STD{i}.vel.data,0,3);
    end
    
    if dim_red_data_yn
        STD{i}.neural_data.dim_red.mn = mean(STD{i}.neural_data.dim_red.scores,3);
        STD{i}.neural_data.dim_red.sd = std(STD{i}.neural_data.dim_red.scores,0,3);
        
        STD{i}.emg_data.dim_red.mn = mean(STD{i}.emg_data.dim_red.scores,3);
        STD{i}.emg_data.dim_red.sd = std(STD{i}.emg_data.dim_red.scores,0,3);
    end
    
    % add bin size
    STD{i}.bin_size         = bin_size;
    % add trial time averaging method
    STD{i}.avg_method       = trial_norm; 
    % add start and end words
    STD{i}.w_i              = w_i;
    STD{i}.w_f              = w_f;
    % create time vector
    STD{i}.t                = 0:bin_size:(round(nbr_bins)-1)*bin_size;
    
    % add info about emgs and neurons
    STD{i}.neural_data.neural_chs = neural_chs;
    STD{i}.emg_data.emg_names = binned_data.emgguide(emg_chs);
end


% -------------------------------
% create a last field in the STD struct that has all the targets
% concatenated
ptr                         = nbr_targets+1;

% ToDo: try to make this code more compact/elegant
aux_fr                      = STD{1}.neural_data.fr;
aux_fr_m                    = STD{1}.neural_data.mn;
aux_fr_sd                   = STD{1}.neural_data.sd;

if isfield(cropped_binned_data,'smoothedspikerate')
    aux_smoothed_fr         = STD{1}.neural_data.smoothed_fr;
    aux_smoothed_fr_m       = STD{1}.neural_data.smoothed_fr_mn;
    aux_smoothed_fr_sd      = STD{1}.neural_data.smoothed_fr_sd;
end

aux_emg                     = STD{1}.emg_data.emg;
aux_emg_m                   = STD{1}.emg_data.mn;
aux_emg_sd                  = STD{1}.emg_data.sd;

if ~isempty(cropped_binned_data.cursorposbin)
    aux_pos                 = STD{1}.pos.data;
    aux_pos_m               = STD{1}.pos.mn;
    aux_pos_sd              = STD{1}.pos.sd;
    aux_vel                 = STD{1}.vel.data;
    aux_vel_m               = STD{1}.vel.mn;
    aux_vel_sd              = STD{1}.vel.sd;
end

if dim_red_data_yn
    aux_neural_scores       = STD{1}.neural_data.dim_red.scores;
    aux_neural_scores_m     = STD{1}.neural_data.dim_red.mn;
    aux_neural_scores_sd    = STD{1}.neural_data.dim_red.sd;
    aux_emg_scores          = STD{1}.emg_data.dim_red.scores;
    aux_emg_scores_m        = STD{1}.emg_data.dim_red.mn;
    aux_emg_scores_sd       = STD{1}.emg_data.dim_red.sd;
end

aux_bin_indx                = STD{1}.bin_indx;

aux_bin_indx_p_trial        = STD{1}.bin_indx_p_trial;


% concatenate trials in 'aux' vars
for i = 2:nbr_targets
    aux_fr                  = cat(3,aux_fr,STD{i}.neural_data.fr);
    aux_fr_m                = cat(1,aux_fr_m,STD{i}.neural_data.mn);
    aux_fr_sd               = cat(1,aux_fr_sd,STD{i}.neural_data.sd);
    
    if isfield(cropped_binned_data,'smoothedspikerate')
        aux_smoothed_fr     = cat(3,aux_smoothed_fr,STD{i}.neural_data.smoothed_fr);
        aux_smoothed_fr_m   = cat(1,aux_smoothed_fr_m,STD{i}.neural_data.smoothed_fr_mn);
        aux_smoothed_fr_sd  = cat(1,aux_smoothed_fr_sd,STD{i}.neural_data.smoothed_fr_sd);
    end
    
    aux_emg                 = cat(3,aux_emg,STD{i}.emg_data.emg);
    aux_emg_m               = cat(1,aux_emg_m,STD{i}.emg_data.mn);
    aux_emg_sd              = cat(1,aux_emg_sd,STD{i}.emg_data.sd);
    
    if ~isempty(cropped_binned_data.cursorposbin)
        aux_pos             = cat(3,aux_pos,STD{i}.pos.data);
        aux_pos_m           = cat(1,aux_pos_m,STD{i}.pos.mn);
        aux_pos_sd          = cat(1,aux_pos_sd,STD{i}.pos.sd);
        aux_vel             = cat(3,aux_vel,STD{i}.vel.data);
        aux_vel_m           = cat(1,aux_vel_m,STD{i}.vel.mn);
        aux_vel_sd          = cat(1,aux_vel_sd,STD{i}.vel.sd);
    end
    
    if dim_red_data_yn
        aux_neural_scores   = cat(3,aux_neural_scores,STD{i}.neural_data.dim_red.scores);
        aux_neural_scores_m = cat(1,aux_neural_scores_m,STD{i}.neural_data.dim_red.mn);
        aux_neural_scores_sd = cat(1,aux_neural_scores_sd,STD{i}.neural_data.dim_red.sd);
        aux_emg_scores      = cat(3,aux_emg_scores,STD{i}.emg_data.dim_red.scores);
        aux_emg_scores_m    = cat(1,aux_emg_scores_m,STD{i}.emg_data.dim_red.mn);
        aux_emg_scores_sd   = cat(1,aux_emg_scores_sd,STD{i}.emg_data.dim_red.sd);
    end
    
    aux_bin_indx            = cat(2,aux_bin_indx,STD{i}.bin_indx);
    aux_bin_indx_p_trial    = cat(2,aux_bin_indx_p_trial,STD{i}.bin_indx_p_trial);
end


% and now add them to the STD struct

STD{ptr}.neural_data.fr     = aux_fr;
STD{ptr}.neural_data.mn     = aux_fr_m;
STD{ptr}.neural_data.sd     = aux_fr_sd;

if isfield(cropped_binned_data,'smoothedspikerate')
    STD{ptr}.neural_data.smoothed_fr    = aux_smoothed_fr;
    STD{ptr}.neural_data.smoothed_fr_mn = aux_smoothed_fr_m;
    STD{ptr}.neural_data.smoothed_fr_sd = aux_smoothed_fr_sd;
end

STD{ptr}.emg_data.emg       = aux_emg;
STD{ptr}.emg_data.mn        = aux_emg_m;
STD{ptr}.emg_data.sd        = aux_emg_sd;

if ~isempty(cropped_binned_data.cursorposbin)
    STD{ptr}.pos.data       = aux_pos;
    STD{ptr}.pos.mn         = aux_pos_m;
    STD{ptr}.pos.sd         = aux_pos_sd;
    STD{ptr}.vel.data       = aux_vel;
    STD{ptr}.vel.mn         = aux_vel_m;
    STD{ptr}.vel.sd         = aux_vel_sd;
end

if dim_red_data_yn
    STD{ptr}.neural_data.dim_red.scores     = aux_neural_scores;
    STD{ptr}.neural_data.dim_red.mn         = aux_neural_scores_m;
    STD{ptr}.neural_data.dim_red.sd         = aux_neural_scores_sd;
    STD{ptr}.emg_data.dim_red.scores        = aux_emg_scores;
    STD{ptr}.emg_data.dim_red.mn            = aux_emg_scores_m;
    STD{ptr}.emg_data.dim_red.sd            = aux_emg_scores_sd;
end


% add bin size
STD{ptr}.bin_size           = bin_size;
% add trial time averaging method
STD{ptr}.avg_method         = trial_norm; 
% add start and end words
STD{ptr}.w_i                = w_i;
STD{ptr}.w_f                = w_f;
% create time vector
STD{ptr}.t                  = 0:bin_size:(size(STD{ptr}.neural_data.mn,1)-1)*bin_size;
STD{ptr}.t_single_trial     = 0:bin_size:(size(STD{ptr}.neural_data.fr,1)-1)*bin_size;

% bin numbers in cropped_binned_data
STD{ptr}.bin_indx           = aux_bin_indx;
STD{ptr}.bin_indx_p_trial   = aux_bin_indx_p_trial;



% ------------------------------------------------------------------------
% return variables

single_trial_data.target    = STD;

