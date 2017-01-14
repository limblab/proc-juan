%
% Build decoder with cross-validation
%
%   mfxval_fits = mfxval_decoder( X, y, varagin )
%
%
% Inputs (opt)      : [default]
%   X               : input data (bins x variables)
%   y               : output data (bins x variables)
%   (lags)          : [10] history lags
%   (fold_size)     : [3000] fold size in bins
%   (dec_offset)    : [true] include offset in decoder model
%   (poly_order)    : [2] order of the polynomial in the static
%                       non-linearity
%
% Outputs:
%   mfxval_fits.R2  : R2 of the model fits
%   mfxval_fits.r   : correlation between input data and predictions
%
%


function mfxval_fits = mfxval_decoder( X, y, varargin )


% Default parameters
params              = struct('lags',        10, ...
                            'fold_size',    3000, ... % nbr of bins
                            'dec_offset',   true, ...
                            'poly_order',   2);

% read input parameters and replace default values where necessary
param_names         = fieldnames(params);
for p = reshape(varargin,2,[])
    if any(strcmp(p{1},param_names))
        params.(p{1})   = p{2};
    else
        error([p{1} ' is not a recognized parameter name']);
    end
end

                        
% Get number of folds
nbr_folds           = floor(size(X,1)/params.fold_size);
% And vector with bin number
bin_indx            = 1:size(X,1);

% preallocage matrices for quality of fit metrics
R2                  = zeros(size(y,2),nbr_folds);
r                   = zeros(size(y,2),nbr_folds);

% -------------------------------------------------------------------------
% Build one decoder for each fold
for f = 0:nbr_folds-1
   
    % 1. compute fold start and end, and get indexes for the train and test
    % data windows
    test_data_start = f*params.fold_size + bin_indx(1);
    test_data_end   = test_data_start + params.fold_size - 1;
    test_data_win   = test_data_start:test_data_end;
    train_data_win  = setdiff(bin_indx,test_data_win);
    
    % 2. split the data into training data and test data
    X_test          = X(test_data_win,:);
    X_train         = X(train_data_win,:);
    y_test          = y(test_data_win,:);
    y_train         = y(train_data_win,:);
    
    % 3. build decoder with test data
    if ~params.dec_offset
        H           = filMIMO3(X_train,y_train,params.lags,1,1);
    else
        H           = filMIMO4(X_train,y_train,params.lags,1,1);
    end
    
    % 4. predict test data with linear filter
    if ~params.dec_offset
        [y_pred,~,y_test_new] = predMIMO3(X_test,H,1,1,y_test);
    else
        [y_pred,~,y_test_new] = predMIMO4(X_test,H,1,1,y_test);
    end
    
    % 5. add nonlinearity if specified
    if params.poly_order > 0
        Hnl         = zeros(params.poly_order+1,size(y,2));
        for v = 1:size(y,2)
            % compute the non-linearity
            Hnl(:,v) = WienerNonlinearity(y_pred(:,v),y_test_new(:,v),...
                        params.poly_order);
            % and predict the data with the cascade decoder
            y_pred(:,v) = polyval(Hnl(:,v),y_pred(:,v));
        end
    end
    
    % 6. get quality of fit metrics
    R2(:,f+1)       = CalculateR2(y_pred,y_test_new);
    r(:,f+1)        = calc_r(y_pred,y_test_new);
end



% return variables
mfxval_fits.R2      = R2;
mfxval_fits.r       = r;