%
% To build Wiener cascade decoders given any set of inputs or outputs. This
% is the function used by the wrappers build_neuron_decoder,
% build_PC_decoder, build_dPC_decoder
%
%   dec_out = build_decoder( X, y, params ) 
%
%
% Inputs:
%   X           : array or cell array with input data 
%   y           : array or cell array with input data
%   params      : decoder parameters --see wrapper functions for details
%
%

function dec_out = build_decoder( X, y, params )


% ------------------------------------
% 0. Train on all the data

% see if we have passed a single or mutliple datasets
if iscell(X)
    nbr_decs        = numel(X);
else
    nbr_decs        = 1;
    % conver input matrix into cell array -for flexibility when computing
    % multiple decoders
    X_old           = X;
    y_old           = y;
    clear X y
    X{1}            = X_old;
    y{1}            = y_old;
end


% ------------------------------------
% 1. Train on all the data

% 1.a. Find the linear model
if ~params.dec_offset
    for b = 1:nbr_decs
        H{b}        = filMIMO3(X{b},y{b},params.lags,1,1);
    end
else
    for b = 1:nbr_decs
        H{b}        = filMIMO4(X{b},y{b},params.lags,1,1);
    end
end

% 1.b. Predict the data with the linear model
if ~params.dec_offset
    for b = 1:nbr_decs
        [y_pred{b},X_all_new{b},y_all_new{b}] = predMIMO3(X{b},H{b},1,1,y{b});
    end
else
    for b = 1:nbr_decs
        [y_pred{b},X_all_new{b},y_all_new{b}] = predMIMO4(X{b},H{b},1,1,y{b});
    end
end

% 1.c. Add non-linearity, if specified
if params.poly_order > 0 
    for b = 1:nbr_decs
        for v = 1:size(y,2)
            % compute the non-linearity
            Hnl{b}(:,v) = WienerNonlinearity(y_pred{b}(:,v),y_all_new{b}(:,v),...
                params.poly_order);
            % and predict the data with the cascade decoder
            y_pred{b}(:,v) = polyval(Hnl{b}(:,v),y_pred{b}(:,v));
        end
    end
end

% 1.d. Get quality of fit
for b = 1:nbr_decs
    R2{b}           = CalculateR2(y_pred{b},y_all_new{b});
    r{b}            = calc_r(y_pred{b},y_all_new{b});
end


% ------------------------------------
% 2. Do multifold cross-validation

if params.xval_yn
    for b = 1:nbr_decs
        mfxval_fits{b} = mfxval_decoder( X{b}, y{b}, 'lags', params.lags, ...
                'fold_size', params.fold_size/params.bin_size, ...
                'dec_offset', params.dec_offset, ...
                'poly_order', params.poly_order );
    end
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% D. OUTPUT DATA

% from training on all
for b = 1:nbr_decs
    dec_out.all{b}.H       = H{b};
    if exist('Hnl','var')
        dec_out.all{b}.Hnl = Hnl{b};
    end
    dec_out.all{b}.r       = r{b};
    dec_out.all{b}.R2      = R2{b};
    if params.return_preds
        dec_out.all{b}.y_pred = y_pred{b};
        dec_out.all{b}.X_ds = X_all_new{b};
        dec_out.all{b}.y   = y_all_new{b};
    end
end

% cross-validated results
for b = 1:nbr_decs
    if params.xval_yn
        dec_out.mfxval{b}.R2 = mfxval_fits{b}.R2;
        dec_out.mfxval{b}.r = mfxval_fits{b}.r;
    end
end

% save parameters
dec_out.params      = params;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% E. PLOTS

if params.plot_yn
    figure,
    marg_leg        = {'task','target','time','task/tgt int'};
    for b = 1:nbr_decs
        if params.xval_yn
            subplot(1,2*nbr_decs,2*b)
            hold on
            errorbar(mean(mfxval_fits{b}.R2,2),std(mfxval_fits{b}.R2,0,2),'marker','none',...
                'color','k','linewidth',2,'linestyle','none')
            bar(mean(mfxval_fits{b}.R2,2))
            ylim([0 1.1]),xlim([0 size(y{b},2)+1])
            ylabel('xval R2'),xlabel('muscle #')
            set(gca,'TickDir','out','FontSize',14)
            box off
            title(marg_leg(params.margs))
            subplot(1,2*nbr_decs,2*b-1)
        end
            bar(R2{b},'FaceColor',[.6 .6 .6])
            ylim([0 1.1]),xlim([0 size(y{b},2)+1])
            ylabel('R2'),xlabel('muscle #')
            set(gca,'TickDir','out','FontSize',14)
            title(marg_leg(params.margs));
            box off
    end
end