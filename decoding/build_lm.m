%
% Build linear models with no history. Optionally, cross-validate the
% predictions. If want to predict several variables, model matrix W is a
% series of MISO models. The model has the form: y = W·X + W_offset
%
%   function [W, W_offset, stats] = build_lm( X, y, xval_yn, fold_length, ...
%                                       xval_plots_yn )
% 
%
% Inputs (opt)          : [defaults]
%   X                   : input data matrix (time-by-inputs)
%   y                   : output data matrix (time-by-outputs)
%   (xval_yn)           : [false] cross-validation flag
%   (fold_length)       : [3000] fold length (nbr of samples)
%   (xval_plots_yn)     : [false] plot xval stuff
%
% Outputs:
%   W                   : linear model
%   W_offset            : offsets for the model
%   stats               : fit stats - adjusted R^2. Without
%                           cross-validation is a single column matrix,
%                           with the R^2 for each output variable in a
%                           different row. With cross-validation it has
%                           nbr_folds + 1 columns, where the first column
%                           is the R^2 using all the data, and the reminder
%                           the R^2 for each variable and each fold
%
%

function [W, W_offset, stats] = build_lm( X, y, varargin )


% -------------------------------------------------------------------------
% read fcn inputs

if nargin == 2
    xval_yn             = false;
    xval_plots_yn       = false;
end
if nargin == 3
    xval_yn             = varargin{1};
    fold_length         = 3000;
    xval_plots_yn       = false;
end
if nargin == 4
    xval_yn             = varargin{1};
    fold_length         = varargin{2};
    xval_plots_yn       = false;
end
if nargin == 5
    xval_yn             = varargin{1};
    fold_length         = varargin{2};
    xval_plots_yn       = varargin{3};
end


% get nbr of output variables we want to predict, and number of input
% variables (predictors)
nbr_outputs             = size(y,2);
nbr_inputs              = size(X,2);


% -------------------------------------------------------------------------
% build the model

% get nbr of folds, if doing multifold cross-validation
if xval_yn
    nbr_folds           = floor( size(X,1) / fold_length );
end

% preallocate model matrix, offset matrix, and matrix to store the adjusted
% R2
W                       = zeros( nbr_outputs, nbr_inputs );
W_offset                = zeros( nbr_outputs, 1 );
if ~xval_yn
    stats               = zeros( nbr_outputs, 1 );
else
    stats               = zeros( nbr_outputs, 1+nbr_folds );
end


% -------------------------------------------------------------------------
% build MISO models using all the data

for v = 1:nbr_outputs

    % fit & store model
    aux_lm              = fitlm( X, y(:,v) );
    W(v,:)              = table2array(aux_lm.Coefficients(2:end,1))';
    % store the intercept (offset)
    W_offset(v)         = table2array(aux_lm.Coefficients(1,1));

    % and the stats
    stats(v,1)          = aux_lm.Rsquared.Adjusted; 
end


% ------------------------------------
% do crossvalidation, if specified

if xval_yn
    
    % preallocate matrix to store all models
    W_xval              = zeros( nbr_outputs, nbr_inputs, nbr_folds );
    W_offset_xval       = zeros( nbr_outputs, 1, nbr_folds );
    
    % do for each fold
    for f = 1:nbr_folds

        % data for testind the model (1 fold)
        start_test      = 1 + (f-1)*fold_length;
        end_test        = f*fold_length;
        bins_test       = start_test:end_test;
        % data for building the model (the rest)
        bins_build      = 1:size(X,1);
        bins_build(ismember(bins_build,bins_test)) = [];
        
        % create vector with data for building and testing
        X_test          = X(bins_test,:);
        y_test          = y(bins_test,:);
        X_build         = X(bins_build,:);
        y_build         = y(bins_build,:);
        
        % build MISO models
        for v = 1:nbr_outputs
           
            % build model with 'build data'
            aux_lm      = fitlm( X_build, y_build(:,v) );
            aux_W       = table2array(aux_lm.Coefficients(2:end,1))';
            aux_W_offset = table2array(aux_lm.Coefficients(1,1));
            
            % predict 'test data' with the model
            y_hat       = repmat(aux_W_offset,1,fold_length) + aux_W*X_test';
            this_R2     = CalculateR2( y_hat', y_test(:,v) );
            
            % save model params
            W_xval(v,:,f)       = aux_W;
            W_offset_xval(v,1,f) = aux_W_offset;
            
            % compute R^2
            stats(v,f+1)  = this_R2; 
        end
    end
end



% -------------------------------------------------------------------------
% PLOT crossvalidation stuff (optional)

if xval_plots_yn
    colors                  = distinguishable_colors(nbr_outputs);
    
    % plot the predictions for each fold
    figure,hold on
    for v = 1:nbr_outputs
        plot(stats(v,2:end),'color',colors(v,:),'marker','.','markersize',24)
        % mean w errorbars
        this_sd             = std(stats(v,2:end));
        this_mn             = mean(stats(v,2:end));
        plot([-1 -1],[this_mn-this_sd,this_mn+this_sd],'linewidth',2,'color',colors(v,:))
        plot(-1,this_mn,'linewidth',2,'color',colors(v,:),'marker','o','markersize',8)
        plot(0,stats(v,1),'v','color',colors(v,:),'markersize',12)
    end
    ylim([0 1]),xlim([-2, nbr_folds+1])
    ylabel('adjusted R^2'),xlabel('fold number')
    set(gca,'TickDir','out','FontSize',14)
    
    % plot the decoders weights for each fold
    nbr_rows                = floor(sqrt(nbr_outputs));
    nbr_cols                = ceil(nbr_outputs/nbr_rows);
    colors2                 = distinguishable_colors(nbr_inputs);
    figure,
    for v = 1:nbr_outputs
        subplot(nbr_rows,nbr_cols,v), hold on
        for i = 1:nbr_inputs
            plot(squeeze(W_xval(v,i,:)),'color',colors2(i,:),'marker','.','markersize',16)
        end
        set(gca,'TickDir','out','FontSize',14)
        ylabel(['output #' num2str(v)])
        
        if v > (nbr_rows-1)*nbr_cols
            xlabel('fold')
        end
        if v == 1, legend('show'), end
    end

    % and the differences in decoder weight as fcn of the fold number
    % - normalized by the abs of the weight
    figure,
    for v = 1:nbr_outputs
        subplot(nbr_rows,nbr_cols,v), hold on
        for i = 1:nbr_inputs
            plot((squeeze(W_xval(v,i,:))-repmat(W_xval(v,i,1),nbr_folds,1))/abs(repmat(W_xval(v,i,1),nbr_folds,1)),...
                'color',colors2(i,:),'marker','.','markersize',16)
        end
        set(gca,'TickDir','out','FontSize',14)
        ylabel(['norm weight change - output #' num2str(v)])        
        if v > (nbr_rows-1)*nbr_cols
            xlabel('fold')
        end
    end
    
    % and the differences in decoder weight as fcn of the fold number
    % - and the real change
    figure,
    for v = 1:nbr_outputs
        subplot(nbr_rows,nbr_cols,v), hold on
        for i = 1:nbr_inputs
            plot((squeeze(W_xval(v,i,:))-repmat(W_xval(v,i,1),nbr_folds,1)),...
                'color',colors2(i,:),'marker','.','markersize',16)
        end
        set(gca,'TickDir','out','FontSize',14)
        ylabel(['weight change - output #' num2str(v)])        
        if v > (nbr_rows-1)*nbr_cols
            xlabel('fold')
        end
    end
end

