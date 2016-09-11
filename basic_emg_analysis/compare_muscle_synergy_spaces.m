%
% Compare muscle spaces, i.e. the hyperplanes defined using PCA or NMF,
% across different tasks. These hyperplanes will be defined as p x m
% matrices, where m is the number of EMGs, and p the dimensionality of the
% muscle hyperplane, i.e. the number of PCs or factors.
%
%   function muscle_space_comp = compare_muscle_spaces( data_struct, varargin ) 
%
%
% Inputs (opt)          : [default]
%   data_struct         : can be a 'single_trial_data' cell array (obtained
%                           with get_single_trial_data), or a 'dataset'
%                           struct (obtained with batch_preprocess_dim_red_data) 
%                           with a field that is a cell array of
%                           dim_red_emg' structs (it can also be any struct
%                           with a cell array with 'dim_red_emg' fields)
%   (nbr_pcs)           : ['empty'] the number of PCs that define the
%                           muscle space. Has to be passed if
%                           dimensionality reduction was done with PCA
%
% Output:
%   muscle_space_comp   : struct with results
%
%

function muscle_space_comp = compare_muscle_synergy_spaces( data_struct, varargin ) 


% -------------------------------------------------------------------------
% see what type of data struct we have and store dim reduction matrices in
% W_cell

% see if it's a dataset type of struct (or a just a cell array of
% dim_red_emg) fields
if isstruct(data_struct)
    if ~isfield(data_struct,'dim_red_emg'), error('data_struct does not have dim_red_emg data'); end
    % get dim red matrices
    W_cell              = cellfun( @(x) x.w, data_struct.dim_red_emg, 'UniformOutput', false );
    dim_red_method      = unique(cellfun( @(x) x.method, data_struct.dim_red_emg, 'UniformOutput', false ));
% or a single_trial_data struct
elseif iscell(data_struct)
    if ~sifield(data_struct{1},'target'), error('data_struct does not have dim_red_emg data'); end
    % get dim red matrices
    W_cell              = cellfun( @(x) x.target{end}.emg_data.dim_red.w, data_struct, ...
                            'UniformOutput', false );
    % and double-check dim red method
    dim_red_method      = unique( cellfun( @(x) x.target{end}.emg_data.dim_red.method, data_struct, ...
                            'UniformOutput', false ) );
end


% -------------------------------------------------------------------------
% read inputs -- if doing PCA, we need the number of dimensions that define
% the muscle synergy space. Otherwise ignore it and give a warning
if nargin == 2
    switch dim_red_method{1}
        case 'pca'
            nbr_syn     = varargin{1};
        case 'nmf'
            disp('the nbr. of synergies (param 2) will be disregarded for NMF')
    end
elseif nargin == 1
    if dim_red_method{1} == 'pca'
        error('PCA was used to get the muscle synergies, so you need to pass the nbr. of synergies as input param 2')
    end
end


% if doing PCA, cut the W_cell matrices to the desired number of dimensions
if dim_red_method{1} == 'pca'
    W_cell          = cellfun( @(x) x(:,1:nbr_syn), W_cell, 'UniformOutput', false );
end

% get nbr of muscle synergies, if using NMF
if ~exist('nbr_syn','var')
    nbr_syn         = size(W_cell{1},2);
end


% -------------------------------------------------------------------------
% compare the muscle synergy spaces

% get nbr of tasks
nbr_bdfs                = length(W_cell);
% and all possible combinations
comb_bdfs               = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs           = size(comb_bdfs,1);


% preallocate matrix for storing results
princ_ang               = zeros(nbr_comb_bdfs, nbr_syn);


% do for all combinations
for i = 1:nbr_comb_bdfs
    princ_ang(i,:)      = principal_angles( W_cell{comb_bdfs(i,1)}, W_cell{comb_bdfs(i,2)} );
end


% return var --struct for compatibility with "neural" code
muscle_space_comp.princ_ang = princ_ang;