%
% Retrieve the latent variables computed with dPCA
%
%
%


function lat_vars = get_lat_vars_dPCA( FR_avg, W, varargin )


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


% -------------------------------------------------------------------------
% reshape to a neuron by time matrix (with all conditions and targets
% concatenated)
FR_conc         = FR_avg(:,:)';
FR_cen          = bsxfun(@minus, FR_conc, mean(FR_conc));
FR_avg_cen      = bsxfun(@minus, FR_avg, mean(FR_conc)');

data_dim        = size(FR_avg);

Z               = FR_cen*W;

% -------------------------------------------------------------------------
% Reshape back to get the latent variables per marginalization 

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