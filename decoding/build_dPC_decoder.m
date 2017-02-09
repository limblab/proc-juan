%
% To build dPC to EMG/force/pos/vel decoders. By default input
% marginalizations are: 1: 'task', 2: 'target', 3: 'time', 4: 'task/target
%
%   dec_out = build_dPC_decoder( ds, dpc, varargin ) 
%
% 
% Inputs (opt)      : [default]
%   ds              : dataset structure, with binned_data and
%                       single_trial_data structures
%   dpc             : dPCA_data structure
%   (w_i)           : ['ot_on'] word that defines the start of the
%                       trial-related data window 
%   (w_f)           : ['R'] word that defines the end of the trial-related
%                       data window 
%   (lags)          : [10] history lags
%   (margs)         : dPCA marginalizations (1: task, 2: target, 3: time,
%                       4: task/target interact.) -- can be multiple
%   (nbr_dpcs_p_marg) : [] number of dPCs per marginalization to use as
%                       decoder inputs. If empty it will use all
%   (dpcs)          : [] what dPCs to use. Incompatible with specifying the
%                       marginalization(s) to use in 'margs'
%   (xval_yn)       : [true] do multifold cross-validation?
%   (fold_size)     : [60] fold size in s
%   (smooth_spikes) : [true] used smoothed or 'raw' spikes
%   (dec_offset)    : [true] include an offset in the decoder? --note that
%                       the spikes are mean-subtracted
%   (poly_order)    : [2] order of the polynomial in the static
%                       non-linearity
%   (output)        : ['emg'] output variable
%   (dec_p_task)    : [false] build a decoder for each task, or one for all
%                       combined tasks
%   (return_preds)  : [true] return input and output data and predictions
%   (plot_yn)       : [false] summary plots
%
% Outputs:
%   dec_out         : structure with decoder weights, predicted data and
%                       quality of fit metrics
%
%
% ToDo:
%   - need to implement force and kinematic predictions
%
%


function dec_out = build_dPC_decoder( ds, dpc, varargin ) 

% Default parameters
params                  = struct('w_i',             'ot_on', ...
                                'w_f',              'R', ...
                                'lags',             10, ...
                                'margs',            [], ...
                                'nbr_dpcs_p_marg',  [], ...
                                'dpcs',             [], ...
                                'xval_yn',          true, ...
                                'fold_size',        60, ...     % in s
                                'smooth_spikes',    true, ...
                                'dec_offset',       true, ...
                                'poly_order',       2, ...
                                'output',           'emg', ...
                                'dec_p_task',       false, ...
                                'return_preds',     true, ...
                                'plot_yn',          false);


% read input parameters and replace default values where necessary
param_names         = fieldnames(params);
for p = reshape(varargin,2,[])
    if any(strcmp(p{1},param_names))
        params.(p{1})   = p{2};
    else
        error([p{1} ' is not a recognized parameter name']);
    end
end


% check that the user hasn't mistakenly passed the dPC numbers and the
% marginalizations to use, as those are not compatible
if ~isempty(params.margs) && ~isempty(params.dpcs)
    error('Cannot choose the marginalization(s) and the dPCs to use as predictors at the same time');
end


% get number of tasks (bdfs)
nbr_bdfs            = length( ds.binned_data );
% get bin size
params.bin_size     = ds.stdata{1}.target{1}.bin_size;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% A. PREPARE THE NEURAL DATA

% 1. Take the binned neural data from the chanels we want 
for b = 1:nbr_bdfs
    % take smoothed data or square root transformed 'raw' spikes 
    if params.smooth_spikes
        sp{b}       = ds.binned_data{b}.smoothedspikerate(:,ds.neural_chs);
    else
        sp{b}       = sqrt(ds.binned_data{b}.spikeratedata(:,ds.neural_chs));
    end
    % subtract the mean FR for each channel
    for n = 1:size(sp{b},2)
        sp{b}(:,n)  = sp{b}(:,n) - ones(size(sp{b},1),1)*mean(sp{b}(:,n),1);
    end
end


% 2. Project them onto the the dPCA axes
for b = 1:nbr_bdfs
    lv{b}           = sp{b}*dpc.W;
end


% 3.a. Choose the marginalizations that we want. By default, the
% marginalizations are 1: 'task', 2: 'target', 3: 'time', 4: 'task/target
% interact'
if ~isempty(params.margs)
    
    % find what dPCs correspond to the input marginalizations
    lvs_marg        = [];
    for i = 1:length(params.margs)
        % see if we want to keep them all
        if isempty(params.nbr_dpcs_p_marg)
            lvs_marg = [lvs_marg , find(dpc.which_marg==params.margs(i))];
        % or just the first few
        else
            lvs_marg = [lvs_marg , find(dpc.which_marg==params.margs(i),...
                        params.nbr_dpcs_p_marg)];
        end
    end
    
    % and keep only those dPCs
    for b = 1:nbr_bdfs
        lv{b}    	= lv{b}(:,lvs_marg);
    end
end

% 3.b. Alternatively, the user may have specified what dPCs to use as
% predictors
if ~isempty(params.dpcs)
    
    for b = 1:nbr_bdfs
        lv{b}     	= lv{b}(:,params.dpcs);
    end
end


% 4. Duplicate and shift the binned data
for b = 1:nbr_bdfs
    lv_ds{b}        = DuplicateAndShift(lv{b},params.lags);
end


% 5. Cut the neural data and the variables to be predicted between the
% desired words 
for b = 1:nbr_bdfs
    % find cutting times
    cut_b{b}        = find_cutting_times(ds.binned_data{b}.trialtable, ...
                                        'task',     ds.labels{b},...
                                        'word_i',   params.w_i,...
                                        'word_f',   params.w_f,...
                                        'bin_size', params.bin_size);
                                    
    % and cut !!!
    lv_ds_cut{b}    = cut_data(lv_ds{b},cut_b{b});
    
    % and store cutting times
    % TODO !!!
end


% 6. Create a matrix with the data from all the trials concatenated
if ~params.dec_p_task
    % all tasks concatenated
    X_all{1}        = cell2mat(lv_ds_cut');
else
    % each task separately
    for b = 1:nbr_bdfs
        X_all{b}    = lv_ds_cut{b};
    end
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% B. PREPARE THE DATA TO BE PREDICTED

% 1. Get the data we want
for b = 1:nbr_bdfs
    switch params.output
        case 'emg'
            y{b}    = ds.binned_data{b}.emgdatabin;
        case 'pos'
            y{b}    = ds.binned_data{b}.cursorposbin;
        case 'force'
            y{b}    = ds.binned_data{b}.forcedatabin;
    end
end


% 2. Cut the data to be predicted between the desired words
for b = 1:nbr_bdfs
    y_cut{b}        = cut_data(y{b},cut_b{b});
end


% 3. Create a matrix with the data from all trials concatenated
if ~params.dec_p_task
    % all tasks concatenated
    y_all{1}        = cell2mat(y_cut');
else
    % each task separately
    for b = 1:nbr_bdfs
        y_all{b}    = y_cut{b};
    end
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% C. BUILD DECODER

dec_out = build_decoder( X_all, y_all, params );


% add input field to params.
dec_out.params.input = 'dPCs';