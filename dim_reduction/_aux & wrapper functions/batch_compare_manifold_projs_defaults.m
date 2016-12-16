%
% Params for 'batch_compare_manifold_projs'
%
%   function params = batch_compare_manifold_projs_defaults( varargin )
%
%
% 'params' is a struct with fields:
%
% Field                 : [default]
%   target              : ['all_conc'] target number
%   dim_manifold        : [20] number of dimensions of the neural manifold
%   P_thr               : [0.001] P value for assessing the significance of
%                           the canonical correlations using Matlab's
%                           built-in test 
%   do_bootstrap        : [true] do bootstrapping to assess the
%                           significance of the canonical correlations     
%   prctile_bootstrap   : [99] the percentile CC that will be chosen as
%                           significance threshold
%   time_win            : [see code] time window that wil be used for the
%                           analysis (t_init, t_end) per session
%   plot_p_session      : [false] whether to plot one figure per session or
%                           not, besides the summary plots
%
%


function params = batch_compare_manifold_projs_defaults( varargin )


params_defaults     = struct( ...
                        'target',           'all_conc', ...
                        'dim_manifold',     20, ...
                        'P_thr',            0.001, ... 
                        'do_bootstrap',     true, ...
                        'nbr_shuffles_bootstrap',10000,...
                        'prctile_bootstrap',99, ...
                        'time_win',         [0.1 0.7; 0.1 0.7; 0.1 0.7;
                                            0.5 1.1; 0.3 0.9; 0.4 1;
                                            0 0.6; 0 0.6; 0 0.6; 
                                            0.3 0.9; 0.2 0.8], ...
                        'plot_p_session',   false );


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

