%
% Compute principal (canonical) angles  between all possible pairs of
% tasks.
%
%   function princ_angles = principal_angles_all_manifolds( dim_red_FR, dims, labels, varargin )
%
% Inputs (opt)      : [default]
%   dim_red_FR      : cell array with dim_red_FR structs that define the
%                       manifold of each task
%   dims            : 1-by-m vector that defines the dimensions (principal
%                       components) for which principal angles will be
%                       computed
%   labels          : cell array with teh labels for each task
%   (angle_orth)    : [angle P<0.001] smallest angle that will be obtained
%                       generating 1D vectors in a space with the
%                       dimensionality of the manifold
%
% Outputs:
%   princ_angles    : struct with fields:
%       angles      : p-by-m matrix with m principal angles between all p
%                       pairs of tasks
%       svdec       : p-dimensional struct with fields U, S, V, which are
%                       the singular value decomposition of the product of
%                       the eigenvector matrices (Qa'*Qb = U*S*V')
%       labels      : p-by-1 cell array with the labels of all pairs of
%                       tasks
%       angle_orth  : angle_orth
%
%


function princ_angles = principal_angles_all_manifolds( dim_red_FR, dims, labels, varargin )


% -------------------------------------------------------------------------
% read inputs
if nargin == 4
    angle_orth      = varargin{1};
end

% if not, passed, compute the orthogonality angle for P < 0.001 
if ~exist('angle_orth','var')
    angle_orth      = rad2deg(acos(3.3/sqrt(length(datasets{1}.dim_red_FR{1}.chs))));
end


% -------------------------------------------------------------------------
% get some meta info
nbr_bdfs            = length(dim_red_FR);

comb_bdfs           = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs       = size(comb_bdfs,1);


% -------------------------------------------------------------------------
% do


for p = 1:nbr_comb_bdfs
    
    % compute principal (canonical) angles
    [angles, U, S, V]   = principal_angles( dim_red_FR{comb_bdfs(p,1)}.w(:,dims), ...
                            dim_red_FR{comb_bdfs(p,2)}.w(:,dims) );
    
    % store results
    princ_angles.angles(p,:) = angles;
    princ_angles.svdec(p).U = U;
    princ_angles.svdec(p).S = S;
    princ_angles.svdec(p).V = V;
    princ_angles.labels{p} = [labels(comb_bdfs(p,1)) labels(comb_bdfs(p,1))];    
end

% add angle orthgonality
princ_angles.angle_orth = angle_orth;