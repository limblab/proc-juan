%
% Params for 'output_potent_null_spaced'
%
%   function params = output_potent_null_spaces_defaults( varargin )
% 
% 
%

function params = output_potent_null_spaces_defaults( varargin )


params_defaults     = struct( ...
                        'input',                    'dim_red', ... % 'spikes'
                        'target',                   'all_conc', ...
                        'dim_neural_manifold',      20, ...
                        'dim_emg_manifold',         4, ...
                        'neural_to_output_delay',   0.05, ... % positive = neural first
                        'output',                   'emg', ... 
                        'trial_averaged',           true, ...
                        'fold_length',              3000, ...         
                        'detailed_plots_yn',        true );

                    
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
