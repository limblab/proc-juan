%
% Generate confidence level for CCA, by shuffling the variables of one set
% over time.
%
%   function signif_boots = bootstrap_canon_corr( U, V, varargin )
%
% Inputs (opt)      : [default]
%   U               : t-by-n matrix with the first set of variables (t is
%                       time)
%   V               : t-by-n matrix with the second set of variables
%   (nbr_shuffles)  : [1,000] nbr of times it will be repeated
%   perc_signif     : [99] percentile that will be given as significance
%                       threshold
%
% Output:
%   signif_boots    : 1-by-n matrix with the CC for each variable that is
%                       obtained by chance
%
%

function signif_boots = bootstrap_canon_corr( U, V, varargin )


% -------------------------------------------------------------------------
% parameters

% read inputs
if nargin >= 3
    nbr_shuffles    = varargin{1};
end
if nargin == 4
    perc_signif     = varargin{2};
end

% if not passed, use defaults 
if ~exist('nbr_shuffles','var')
    nbr_shuffles    = 1000; 
end
if ~exist('perc_signif','var')
    perc_signif     = 99;
end


% -------------------------------------------------------------------------
% do !

% get number of neural projections
nbr_projs           = size(U,2);

% preallocate matrices
U_shuffled          = zeros(size(U,1),size(U,2));
cc                  = zeros(nbr_projs,nbr_shuffles);


% shuffle the projections (each independently) of one of them and do CCA
for i = 1:nbr_shuffles
    for p = 1:nbr_projs
        U_shuffled(:,p) = U(randperm(size(U,1)),p);
    end
    
    % CCA
    [~,~,cc(:,i)]   = canoncorr(U_shuffled,V);
end


% get significance value per CC
signif_boots        = prctile(cc',perc_signif);

