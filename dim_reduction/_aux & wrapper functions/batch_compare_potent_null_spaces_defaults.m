%
% Params for 'batch_compare_potent_null_spaces_defaults'
%
%   function params = batch_compare_potent_null_spaces_defaults( varargin )
%
% Inputs (opt)      : [default]
%   (params)        : can be a struct of type
%                       batch_compare_potent_null_spaces_defaults or of
%                       type batch_compare_manifold_projs_defaults, in
%                       which case it will be used to check the consistency
%                       of some params like 'dim_neural_manifold', or
%                       'time_win.' It can be also both structs, the order
%                       doesn't matter
%
% Output:
%   params          : batch_compare_potent_null_spaces_defaults struct
%
%


function params = batch_compare_potent_null_spaces_defaults( varargin )


params_defaults     = struct( ...
                        'input',                    'dim_red', ... % 'spikes'
                        'target',                   'all_conc', ...
                        'dim_neural_manifold',      20, ...
                        'dim_emg_manifold',         4, ... % only for PCA !!! --for NMF will use all the existing ones
                        'neural_to_output_delay',   0.05, ... % positive = neural first
                        'output',                   'emg', ... % pos, vel, force, dim_red_emg
                        'trial_averaged',           true, ...
                        'P_thr',                    0.001, ...                         
                        'time_win',                 [0.2 1.2; 0.2 1.2; 0.2 1.2;
                                                    0 0.74; 0 0.74; 0 0.74; 
                                                    0.1 1; 0.1 1] );

                    
% read input params, if passed
if nargin >= 1
    for i = 1:nargin
        if isfield(varargin{i},'dim_manifold')
            proj_params             = varargin{i};
        elseif isfield(varargin{i},'dim_neural_manifold')
            params                  = varargin{i};
        end
    end
else
    params                          = [];
end


% check consistency of params and proj_params; if they are different update
% params to the values in proj_params and give a warning
if exist('proj_params','var')
    if ~isempty(params)
        if params.target ~= proj_params.target
            params.target           = proj_params.target;
            warning(['params.target updated to ' proj_params.target]);
        end
        if params.dim_neural_manifold ~= proj_params.dim_manifold
            params.dim_neural_manifold  = proj_params.dim_manifold; 
            warning(['params.dim_neural_manifold updated to ' num2str(proj_params.dim_manifold)]);
        end
        if sum(sum(params.time_win ~= proj_params.time_win))
            params.time_win         = proj_params.time_win;
            warning('params.tim_win updated to ');
            disp(num2str(proj_params.time_win));
        end
    end
end


% fill default options missing from input argument
all_params_names    = fieldnames(params_defaults);
for i = 1:numel(all_params_names)
    if ~isfield(params,all_params_names(i))
        params.(all_params_names{i}) = params_defaults.(all_params_names{i});
    end
end
