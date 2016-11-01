%
% Compare latent variables, obtained with PCA, after dropping a certain
% percentage of channels, for all tasks in a dataset array. This function
% is a wrapper of 'comp_latent_vars_dropping_ch.m'
%
%   CC_array = batch_comp_latent_vars_dropping_ch( datasets, perc_drop, varargin )
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
%
% Outputs:
%   CCs             : canonical correlations (with size: iteration x
%                       dimension x perc chs)
%
%

function CCs = batch_comp_latent_vars_dropping_ch( datasets, mani_dim, perc_drop, varargin )


% check that the parallel pool is running, otherwise start it
gcp;

% -------------------------------------------------------------------------
% read inputs

if nargin >= 4
    nbr_reps        = varargin{1};
else
    nbr_reps        = 100;
end
    
if nargin == 5        
    pick_chs        = varargin{2};
else
    pick_chs        = true;
end

nbr_datasets        = length(datasets);


% -------------------------------------------------------------------------
% do
  
% if datasets is a cell array
if nbr_datasets > 1

    for d = 1:nbr_datasets
        tasks_this      = length(datasets{d}.binned_data);
        
        for t = 1:length(tasks_this)    
            if pick_chs
                CCs{d}.corr{t} = comp_latent_vars_dropping_ch( datasets{d}.binned_data{t}.smoothedspikerate, ...
                    mani_dim, perc_drop, nbr_reps, datasets{d}.neural_chs );
            else
                CCs{d}.corr{t} = comp_latent_vars_dropping_ch( datasets{d}.binned_data{t}.smoothedspikerate, ...
                    mani_dim, perc_drop, nbr_reps );
            end
        end
        % store some info
        CCs{d}.labels   = datasets{d}.labels;
        CCs{d}.monkey   = datasets{d}.monkey;
        CCs{d}.data     = datasets{d}.monkey;
    end    
% or if there's only one
else
    tasks_this          = lenght(datasets{d}.binned_data);
    
    for t = 1:length(tasks_this)
        if pick_chs
            CCs.corr{t} = comp_latent_vars_dropping_ch( datasets.binned_data{t}.smoothedspikerate, ...
                mani_dim, perc_drop, nbr_reps, datasets.neural_chs );
        else
            CCs.corr{t} = comp_latent_vars_dropping_ch( datasets.binned_data{t}.smoothedspikerate, ...
                mani_dim, perc_drop, nbr_reps );
        end
    end
    % store some info
    CCs.labels      = datasets.labels;
    CCs.monkey      = datasets.monkey;
    CCs.data        = datasets.monkey;
end


% -------------------------------------------------------------------------
% summary plot

if nbr_datasets > 1

else
    
end
