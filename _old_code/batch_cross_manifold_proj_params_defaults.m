%
% Params for 'batch_cross_manifold_proj'
%
%   function params = batch_cross_manifold_proj_params_defaults( varargin )
%

function params = batch_cross_manifold_proj_params_defaults( varargin )


params_defaults     = struct( ...
                        'dim_manifold',     20, ...
                        'plot_p_pair',      false);
                    

% read input params, if passed
if nargin
    params          = varargin{1};
else
    params          = [];
end

% fill default params missing from input argument
all_params_names    = fieldnames(params_defaults);
for i = 1:numel(all_params_names)
    if ~isfield(params,all_params_names(i))
        params.(all_params_names{i}) = params_defaults.(all_params_names{i});
    end
end