%
% To build PC to EMG/force/pos/vel decoders
%
%   dec_out = build_PC_decoder( ds, varargin ) 
%
% 
% Inputs (opt)      : [default]
%   ds              : dataset structure, with binned_data and
%                       single_trial_data structures
%   (w_i)           : ['ot_on'] word that defines the start of the
%                       trial-related data window 
%   (w_f)           : ['R'] word that defines the end of the trial-related
%                       data window 
%   (lags)          : [10] history lags
%   (nbr_pcs)       : [] number of dPCs per marginalization to use as
%                       decoder inputs. If empty it will use all
%   (pcs)          : [] what dPCs to use. Incompatible with specifying the
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


function dec_out = build_PC_decoder( ds, varargin ) 

% Default parameters
params                  = struct('w_i',             'ot_on', ...
                                'w_f',              'R', ...
                                'lags',             10, ...
                                'nbr_pcs',          [], ...
                                'pcs',              [], ...
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
% number of PCs to use, as those are not compatible
if ~isempty(params.pcs) && ~isempty(params.nbr_pcs)
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


% 2. Project them onto the the PCA axes
for b = 1:nbr_bdfs
    lv{b}           = sp{b}*ds.dim_red_FR{b}.w;
end


% 3.a. Choose the PCs that we want --the user can either pass the PC
% numbers or how many
if ~isempty(params.pcs)
    % keep only these PCs
    for b = 1:nbr_bdfs
        lv{b}    	= lv{b}(:,params.pcs);
    end
end


% 3.b. Alternatively, the user may have specified what PCs to use as
% predictors
if ~isempty(params.nbr_pcs)
    % keep only the leading X PCs
    for b = 1:nbr_bdfs
        lv{b}    	= lv{b}(:,1:params.nbr_pcs);
    end
end


% 4. Duplicate and shift the binned data
for b = 1:nbr_bdfs
    lv_ds{b}        = DuplicateAndShift(lv{b},params.lags);
end


% 5. Cut the neural data and the variable to be predicted between the
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
end


% 6. Create a matrix with the data from all the trials concatenated
X_all               = cell2mat(lv_ds_cut');



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


% 3. Create a matrix with the data from all traisl concatenated
y_all               = cell2mat(y_cut');



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% C. BUILD DECODER

dec_out = build_decoder( X_all, y_all, params );


% add input field to params.
dec_out.params.input = 'PCs';