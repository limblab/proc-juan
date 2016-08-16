

function results = batch_cross_manifold_proj( varargin )


% -------------------------------------------------------------------------
% some options that could be turn into fcn params
dim_manifold            = 20; % the neural manifolds are defined by the first 20 dimensions


% -------------------------------------------------------------------------
% read input params

if nargin >=1
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end

if nargin == 2
    angle_results       = varargin{2};
end

clear varargin;


% -------------------------------------------------------------------------
% load data

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end


% -------------------------------------------------------------------------
% get info about monkeys, tasks, and who did what

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% assess how similat the within-manifold and across-manifold projections
% are

for i = 1:meta_info.nbr_monkeys
     
    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
        
        % what dataset are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);

        pc_projs        = transform_and_compare_dim_red_data_all_tasks( ...
                            datasets{dtst}.dim_red_FR, ...
                            datasets{dtst}.smoothed_FR, ...
                            datasets{dtst}.labels, ...
                            datasets{dtst}.neural_chs, ...
                            dim_manifold, ...
                            'min_angle', ...
                            angle_results.data{dtst}.angles.pair_min_angle );
                        
        % store all results
        proj_comp{dtst} = pc_projs; 
    end
end

