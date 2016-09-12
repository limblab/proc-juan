%
% Do PCA/NMF on the concatenated EMGs of many trials.
%
% function dim_red_emg = dim_reduction_muscles( binned_data, varargin )
%
% Inputs (opt)          : [default]
%   data_set            : data_set struct generated with
%                           batch_preprocess_dim_red_data
%   (method)            : ['pca'] dim reduction method ('pca','nmf','none')

%   (labels)            : name of the task in each binned_data struct
%   (nbr_factors)       : nbr of factors we want to extract (only with NMF)
%   (plot_yn)           : [false] summary plot
%
% Outputs
%   dim_red_emg         : cell with fields:
%       w               : 'pca': eigenvectors in columns, 'nmf': weights in
%                           columns (muscles x nbr_factors)
%       eigen           : eigenvalues of w
%       scores          : result of applying w to smoothed_FR
%       t               : time axis for scores
%       chs             : EMG signals included in the analysis
%       method          : method used (stores input)
%
%


function dim_red_emg = dim_reduction_muscles_pooling_tasks( data_set, varargin )


% -------------------------------------------------------------------------
% get inputs --inhereted from dim_reduction_muscles
if nargin == 1
    method              = 'pca';
    chosen_emgs         = 'all';
elseif nargin >= 2
    method              = varargin{1};
end
if nargin >= 3
    labels              = varargin{2};
end
if nargin >= 4
    nbr_factors         = varargin{3};
end
if nargin == 5
    plot_yn             = varargin{4};
end

if ~exist('plot_yn','var')
    plot_yn             = false;
end
if ~exist('labels','var')
    labels              = [];
end
if strcmp(method,'nmf')
    if ~exist('nbr_factors','var')
        error('need to specify nbr of factors for NMF');
    end
else
    nbr_factors         = [];
end

% the EMG channels we want to use is already defined in data_set
chosen_emgs             = data_set.chosen_emgs;

% -------------------------------------------------------------------------
% nbr of tasks in dataset
nbr_bdfs                    = length(data_set.stdata);

% concatenate all the EMGs
conc_emg                    = data_set.stdata{1}.target{end}.emg_data.conc_emg;
for i = 2:nbr_bdfs
    conc_emg                = cat(1,conc_emg,data_set.stdata{i}.target{end}.emg_data.conc_emg);
end

% do !
dim_red_emg                 = dim_reduction_muscles( conc_emg, method, chosen_emgs, ...
                                labels, nbr_factors, plot_yn );
                            
% add time if it doesn't exist
if ~isfield(dim_red_emg,'t')
    aux_t                   = 0:data_set.stdata{1}.target{1}.bin_size:...
                                data_set.stdata{1}.target{1}.bin_size*(size(dim_red_emg.scores,1)-1);
    dim_red_emg.t           = aux_t';
end