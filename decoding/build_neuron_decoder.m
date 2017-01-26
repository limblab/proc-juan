%
% To build neurons to EMG/force/pos/vel decoders, training and validating
% on multiple tasks
%
%   dec_out = build_neuron_decoder( ds, dpc, varargin ) 
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
%   (choose_chs)    : [true] use only the subset of channels specified in
%                       the 'neural_chs' field of 'ds'
%   (xval_yn)       : [true] do multifold cross-validation?
%   (fold_size)     : [60] fold size in s
%   (smooth_spikes) : [true] used smoothed or 'raw' spikes --note that
%                       raw spikes will be square root transformed anyway!!
%   (dec_offset)    : [true] include an offset in the decoder? --note that
%                       the spikes are mean-subtracted
%   (poly_order)    : [2] order of the polynomial in the static
%                       non-linearity
%   (output)        : ['emg'] output variable
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


function dec_out = build_neuron_decoder( ds, varargin ) 

% Default parameters
params                  = struct('w_i',             'ot_on', ...
                                'w_f',              'R', ...
                                'lags',             10, ...
                                'choose_chs',       true, ...
                                'xval_yn',          true, ...
                                'fold_size',        60, ...     % in s
                                'smooth_spikes',    true, ...
                                'dec_offset',       true, ...
                                'poly_order',       2, ...
                                'output',           'emg', ...
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


% get number of tasks (bdfs)
nbr_bdfs            = length( ds.binned_data );
% get bin size
bin_size            = ds.stdata{1}.target{1}.bin_size;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% A. PREPARE THE NEURAL DATA

% 1. Take the binned neural data from the chanels we want 
for b = 1:nbr_bdfs
    % take smoothed data or square root transformed 'raw' spikes 
    if params.smooth_spikes
        sp{b}       = ds.binned_data{b}.smoothedspikerate;
    else
        sp{b}       = sqrt(ds.binned_data{b}.spikeratedata);
    end
    % subtract the mean FR for each channel
    for n = 1:size(sp{b},2)
        sp{b}(:,n)  = sp{b}(:,n) - ones(size(sp{b},1),1)*mean(sp{b}(:,n),1);
    end
end


% 2. Choose the channels that we want -- the subset in 'neural_chs' or all
if params.choose_chs
    for b = 1:nbr_bdfs
        sp{b}       = sp{b}(:,ds.neural_chs);
    end
end


% 3. Duplicate and shift the binned data
for b = 1:nbr_bdfs
    sp_ds{b}        = DuplicateAndShift(sp{b},params.lags);
end


% 4. Cut the neural data and the variable to be predicted between the
% desired words 
for b = 1:nbr_bdfs
    % find cutting times
    cut_b{b}        = find_cutting_times(ds.binned_data{b}.trialtable, ...
                                        'task',     ds.labels{b},...
                                        'word_i',   params.w_i,...
                                        'word_f',   params.w_f,...
                                        'bin_size', bin_size);
                                    
    % and cut !!!
    sp_ds_cut{b}    = cut_data(sp_ds{b},cut_b{b});
end


% 5. Create a matrix with the data from all the trials concatenated
X_all               = cell2mat(sp_ds_cut');



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

% ------------------------------------
% 1. Train on all the data

% 1.a. Find the linear model
if ~params.dec_offset
    H               = filMIMO3(X_all,y_all,params.lags,1,1);
else
    H               = filMIMO4(X_all,y_all,params.lags,1,1);
end

% 1.b. Predict the data with the linear model
if ~params.dec_offset
    [y_pred,X_all_new,y_all_new] = predMIMO3(X_all,H,1,1,y_all);
else
    [y_pred,X_all_new,y_all_new] = predMIMO4(X_all,H,1,1,y_all);
end

% 1.c. Add non-linearity, if specified
if params.poly_order > 0 
    for v = 1:size(y_all,2)
        % compute the non-linearity
        Hnl(:,v)    = WienerNonlinearity(y_pred(:,v),y_all_new(:,v),...
                        params.poly_order);
        % and predict the data with the cascade decoder
        y_pred(:,v) = polyval(Hnl(:,v),y_pred(:,v));
    end            
end

% 1.d. Get quality of fit
R2                  = CalculateR2(y_pred,y_all_new);
r                   = calc_r(y_pred,y_all_new);


% ------------------------------------
% 2. Do multifold cross-validation

if params.xval_yn
    mfxval_fits = mfxval_decoder( X_all, y_all, 'lags', params.lags, ...
                                    'fold_size', params.fold_size/bin_size, ...
                                    'dec_offset', params.dec_offset, ...
                                    'poly_order', params.poly_order );
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% D. OUTPUT DATA

% from training on all
dec_out.all.H       = H;
if exist('Hnl','var')
    dec_out.all.Hnl = Hnl;
end
dec_out.all.r       = r;
dec_out.all.R2      = R2;
if params.return_preds
    dec_out.all.y_pred = y_pred;
    dec_out.all.X_ds = X_all_new;
    dec_out.all.y   = y_all_new;
end

% cross-validated results
if params.xval_yn
    dec_out.mfxval.R2 = mfxval_fits.R2;
    dec_out.mfxval.r = mfxval_fits.r;
end

% save parameters
dec_out.params      = params;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% E. PLOTS

if params.plot_yn
   
    marg_leg        = {'task','target','time','task/tgt int'};
    figure,
    if params.xval_yn
        subplot(122)
        hold on
        errorbar(mean(mfxval_fits.R2,2),std(mfxval_fits.R2,0,2),'marker','none',...
            'color','k','linewidth',2,'linestyle','none')
        bar(mean(mfxval_fits.R2,2))
        ylim([0 1.1]),xlim([0 size(y_all,2)+1])
        ylabel('xval R2'),xlabel('muscle #')
        set(gca,'TickDir','out','FontSize',14)
        box off
        title(marg_leg(params.margs))
        subplot(121)
    end
        bar(R2,'FaceColor',[.6 .6 .6])
        ylim([0 1.1]),xlim([0 size(y_all,2)+1])
        ylabel('R2'),xlabel('muscle #')
        set(gca,'TickDir','out','FontSize',14)
        title(marg_leg(params.margs));
        box off
end