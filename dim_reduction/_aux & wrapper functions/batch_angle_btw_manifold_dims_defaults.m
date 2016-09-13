%
% Params for 'batch_angle_btw_manifold_dims'
%
%   function params = batch_angle_btw_manifold_dims_defaults( varargin )
%
%
% 'params' is a struct with fields:
%
% Field                 : [default]
%   dim_manifold        : [20] number of dimensions of the neural manifold
%   P_thr               : [0.001] P value for assessing the significance of
%                           the principal angles (obtained from the
%                           empirically generated principal angles)
%   empir_angle_dist_file : [see code] absolute path to the file with the
%                           empirical principal angle data
%   nbr_planes_bootstrap : [10000] if it is necessary to generate a new
%                           empirical distribution of principal angles,
%                           this is the number of samples that will be used
%   plot_p_session      : [false] whether to plot one figure per session or
%                           not, besides the summary plots
%
%


function params = batch_angle_btw_manifold_dims_defaults( varargin )


params_defaults     = struct( ...
                        'dim_manifold',         20, ...
                        'P_thr',                0.001, ...
                        'empir_angle_dist_file', ['/Users/juangallego/Documents/NeuroPlast/Data/' ...
                                                    '_Dimensionality reduction/_control analyses/' ...
                                                    'empirical principal angle distributions all datasets.mat'], ...
                        'nbr_planes_bootstrap', 10000, ...
                        'plot_p_session',       false );
                    
                    
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