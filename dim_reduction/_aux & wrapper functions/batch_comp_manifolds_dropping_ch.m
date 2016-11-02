%
% Compare neural manifolds, obtained with PCA, taking pairs subsets with a
% certain % of all the channels, for all tasks in a dataset array. This
% function is a wrapper of 'comp_latent_vars_dropping_ch.m'
%
%   PA_array = batch_comp_manifolds_dropping_ch( datasets, perc_drop, varargin )
%
%
% Inputs (opt)      : [default]
%   datasets        : dataset array
%   mani_dim        : manifold dimensionality (= nbr. latent variables)
%   perc_drop       : percentage of channels to drop (0-1). It can be an
%                       array, in which case the function will return a 3D
%                       matrix
%   (nbr_reps)      : [100] number of combinations
%   (pick_chs)      : [true] only use the neural channels in
%                       datasets{i}.neural_chs
%   (normal)        : ['sqrt'] normalization of the spike firings
%   (ead)           : empirical angle distribution structure
%
% Outputs:
%   PAs             : principal angles (with size: iteration x
%                       dimension x perc chs)
%
%


function PAs = batch_comp_manifolds_dropping_ch( datasets, mani_dim, perc_drop, varargin )


% check that the parallel pool is running, otherwise start it
gcp;

% -------------------------------------------------------------------------
% read inputs

if nargin >= 4
    nbr_reps        = varargin{1};
else
    nbr_reps        = 100;
end
    
if nargin >= 5        
    pick_chs        = varargin{2};
else
    pick_chs        = true;
end

if nargin >= 6
    normal          = varargin{3};
else 
    normal          = 'sqrt';
end

if nargin == 7
    ead             = varargin{4};
else
    % load the empirical angle distribution (significance level)
    ead             = load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets.mat');
end

nbr_datasets        = length(datasets);


% -------------------------------------------------------------------------
% do
  
% if datasets is a cell array
if nbr_datasets > 1

    for d = 1:nbr_datasets
        tasks_this      = length(datasets{d}.binned_data);
        
        for t = 1:tasks_this    
            if pick_chs
                PAs{d}.ang{t}   = comp_manifolds_dropping_ch( datasets{d}.binned_data{t}.smoothedspikerate, ...
                    mani_dim, perc_drop, nbr_reps, datasets{d}.neural_chs, normal, ead );
            else
                PAs{d}.ang{t}   = comp_manifolds_dropping_ch( datasets{d}.binned_data{t}.smoothedspikerate, ...
                    mani_dim, perc_drop, nbr_reps, [], normal, ead );
            end
        end
        % store some info
        PAs{d}.labels   = datasets{d}.labels;
        PAs{d}.monkey   = datasets{d}.monkey;
        PAs{d}.data     = datasets{d}.monkey;
    end    
% or if there's only one
else
    tasks_this          = lenght(datasets{d}.binned_data);
    
    for t = 1:length(tasks_this)
        if pick_chs
            PAs.corr{t} = comp_manifolds_dropping_ch( datasets.binned_data{t}.smoothedspikerate, ...
                mani_dim, perc_drop, nbr_reps, datasets.neural_chs, normal, ead );
        else
            PAs.corr{t} = comp_manifolds_dropping_ch( datasets.binned_data{t}.smoothedspikerate, ...
                mani_dim, perc_drop, nbr_reps, [], normal, ead );
        end
    end
    % store some info
    PAs.labels      = datasets.labels;
    PAs.monkey      = datasets.monkey;
    PAs.data        = datasets.monkey;
end


% -------------------------------------------------------------------------
% summary plot

if nbr_datasets > 1

else
    
end