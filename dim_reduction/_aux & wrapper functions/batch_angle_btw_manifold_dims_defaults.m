%
% Params for 'batch_angle_btw_manifold_dims'
%
%   function params = batch_angle_btw_manifold_dims_defaults( varargin )
%

function params = batch_angle_btw_manifold_dims_defaults( varargin )


params_defaults     = struct( ...
                        'dim_manifold',         20, ...
                        'P_thr',                0.001, ...
                        'P_empir_angle_dist_file', ['/Users/juangallego/Documents/NeuroPlast/Data/' ...
                                                    '_Dimensionality reduction/_control analyses/' ...
                                                    'empirical principal angle distributions all datasets.mat'], ...
                        'nbr_planes_bootstrap', 10000 );
                    
                    
% read input params, if passed
if nargin
    params          = varargin{1};
else
    params          = [];
end
   
% fill default options missing from input argument
all_params_names    = fieldnames(params_defaults);
for i = 1:numel(all_params_names)
    if ~isfield(params,all_params_names(i))
        params.(all_params_names{i}) = params_defaults.(all_params_names{i});
    end
end