%
% Params for 'batch_preprocess_dim_red_data'
%
%   function params = batch_preprocess_dim_red_data_params_defaults( varargin )
%

function params = batch_preprocess_dim_red_data_params_defaults( varargin )


params_defaults     = struct( ...
                        'bin_size',         0.02, ...
                        'kernel_SD',        0.05, ...
                        'transform',        'sqrt', ...
                        'dim_red_method',   'pca', ...
                        'data_pca',         'trial-related', ...% 'trial-related', 'all'
                        'w_i',              'ot_on', ...
                        'w_f',              'R', ...
                        'norm_trial_data',  'min_dur', ... % 'min_dur', 'stretch'
                        'dim_red_emg',      'nmf', ... % 'none', 'pca','nmf'
                        'emg_factors',       4 ); % only used with NMF

                
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