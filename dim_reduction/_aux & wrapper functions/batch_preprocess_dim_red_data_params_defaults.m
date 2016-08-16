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
                        'dim_red_method',   'pca' );

                
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