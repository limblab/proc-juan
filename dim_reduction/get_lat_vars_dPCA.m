%
% Retrieve the latent variables computed with dPCA
%
%   [lat_vars, varargout] = get_lat_vars_dPCA( FR_avg, FR, W, varargin )
%
%
% Inputs (opt)          : [default]
%   FR_avg              : 4-D matrix with trial-averaged neural data
%                           (neurons x task x target x time)
%   FR                  : 5-D matrix with neural data
%                           (neurons x task x target x time x trial)
%   W                   : dPCA weight matrix
%   (options)           : options -see code for names and default values
%
%


function [lat_vars, varargout] = get_lat_vars_dPCA( FR_avg, FR, W, varargin )


% default input parameters
options         = struct('time',   [], ...   
                 'whichMarg',      [], ...
                 'timeEvents',     [], ...
                 'ylims',          [], ...
                 'componentsSignif', [], ...
                 'timeMarginalization', [], ...
                 'legendSubplot',  [], ...
                 'marginalizationNames', [], ...
                 'marginalizationColours', [], ...
                 'explainedVar',   [], ...
                 'numCompToShow',  15, ...
                 'X_extra',        []);

             
% read input parameters and fill non-specified fields with default values
option_names    = fieldnames(options);
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
    end
end
             

% check that don't want to retrieve more components that there are
num_comps       = min(options.numCompToShow, size(W,2));



% =========================================================================
% =========================================================================
% FOR THE TRIAL-AVERAGED DATA

% -------------------------------------------------------------------------
% 1) reshape to a neuron by time matrix (with all conditions and targets
% concatenated)
FR_conc         = FR_avg(:,:)';
FR_cen          = bsxfun(@minus, FR_conc, mean(FR_conc));
FR_avg_cen      = bsxfun(@minus, FR_avg, mean(FR_conc)');

data_dim        = size(FR_avg);

Z               = FR_cen*W;

% -------------------------------------------------------------------------
% 2) Reshape back to get the latent variables per marginalization 

% get what latent variables we want to retrieve
if isempty(options.componentsSignif)
    comps_to_get = 1:num_comps;
else
    comps_to_get = options.componentsSignif; % ToDo: double check this
end

lat_vars        = reshape(Z(:,comps_to_get)', [length(comps_to_get) data_dim(2:end)]);

% Check that there's no X_extra
if ~isempty(options.X_extra)
    error('options.X_extra not implemented yet')
end



% =========================================================================
% =========================================================================
% FOR SINGLE TRIAL DATA (IF PASSED)

if ~isempty(FR)
   
    FR_c_st     = FR(:,:)';
    FR_cen_st   = bsxfun(@minus, FR_c_st, nanmean(FR_c_st));
    
    data_dim_st = size(FR);
    
    Z_st        = FR_cen_st*W;
    Zfull_st    = reshape(Z_st', [size(Z_st,2) data_dim_st(2:end)]);
end


% -------------------------------------------------------------------------
% Return single trial latent variables if there are two outputs

if nargout == 2
    varargout{1} = Zfull_st;
end