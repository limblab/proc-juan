%
% A wrapper of comp_neural_spaces_fcn_dim that first looks for the
% eigenvectors in manifold 2 that are closest (define the smallest angle)
% which each dimension of manifold 1, for dimensions 1:dims_hyper_in_orig,
% before it computes the angle. 
% By default, the function will also look for the eigenvectors in manifold
% 1 that are closest to each eigenvector in manifold 2, for the same
% dimensions, and give XXXXXXXXXXXXXXXXX. This is to avoid ranking according to one arbitrarily selected
% task
%
%
%   function [angles, dim_red_FR, smoothed_FR empir_angle_dist] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( ...
%           bdf, neural_chs, dims_hyper_in_orig, labels, method, varargin ) 
%
%
% Inputs (opt)          : [default]
%   bdf                 : struct with BDFs. Function will calculate the
%                           angle between the hyperplanes from all possible
%                           pairs of tasks
%   neural_chs          : neural channels to be used for the analysis.
%                           Relevant if the fcn has to compute the
%                           hyperplane, which is done using dim_reduction.m
%   dims_hyper_in_orig  : dimensions in hyperplane 1 that will be used to
%                           define the hyperplanes (function will compare
%                           hyperplanes with dimensionality
%                           1:dims_hyper_in_orig 
%   labels              : cell array that defines the labels for each task
%   method              : dimensionality reduction method ('pca','fa').
%                           Only used if the function has to compute them
%   (smoothed_FR)       : [] cell array of smoothed firing rates (computed 
%                           using gaussian_smoothing2.m)                            
%   (dim_red_FR)        : [] cell array of dimensionally-reduced FRs
%                           (computed using dim_reduction.m)
%   (bin_width)         : [] bin size. Only used if user doesn't pass
%                           smoothed_FR and the function has to compute
%                           them
%   (gauss_SD)          : [] SD of the Gaussian kernel used for smoothing.
%                           Only used if user doesn't pass smoothed_FR and
%                           the function has to compute them
%   (transform)         : [] 'sqrt' or 'none': transformation applied on
%                           the binned FRs. Only used if user doesn't pass
%                           smoothed_FR and the function has to compute
%                           them 
%   (resort_eigenv)     : [false] bool that tells the function to resort
%                           all pairs of eigenvectors not by the variance
%                           of the first task, but based on the angle
%                           between pairs (from smallest to largest) 
%   (last_dim)          : last dimension for resort_eigenv
%
%
% Outputs (opt):
%   angles              : cell array with angles between hyperplanes of
%                           dimensionality 1:dims_hyper_in_orig for all
%                           pairs of tasks
%   (dim_red_FR)        : dim_red_FR for each BDF, if the function has to
%                           compute them
%   (smoothed_FR)       : smoothed_FR, if the function has to compute them
%   (empir_angle_dist)  : empirical distributions of angles between
%                           hyperplanes, used to assess significance
%
%
% Usage:
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method ) 
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method, smoothed_FR)
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method, smoothed_FR, dim_red_FR )
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method, smoothed_FR, dim_red_FR,...
%           empir_angle_dist )
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method, bin_width, gauss_SD, ...
%           transform )  
%   function [angles, dim_red_FR, smoothed_FR, empir_angle_dist ] = ...
%           comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
%           dims_hyper_in_orig, labels, method, smoothed_FR, dim_red_FR,
%           resort_eigenv, last_dim )  
%
%
% Notes/ToDo's:
%   - make 'method' an optional argument
%   - code needs to be cleaned out, e.g. get rid of the need to pass s BDF
%



function [angles, dim_red_FR, smoothed_FR, empir_angle_dist] = ...
            comp_neural_spaces_fcn_dim_finding_closest( bdf, neural_chs,...
            dims_hyper_in_orig, labels, method, varargin ) 
                            

        
% -------------------------------------------------------------------------
% path and file with empirical distributions of angles, for assessing
% orthogonality 
%
% ~-~-> ToDo: turn into a parameter


empir_angle_dist_file       = '/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical angle distribution all datasets.mat';

        
% -------------------------------------------------------------------------
% read input parameters


if nargin == 6
    smoothed_FR             = varargin{1};
elseif nargin == 7
    smoothed_FR             = varargin{1};
    dim_red_FR              = varargin{2};
elseif nargin == 8
    if isscalar(varargin{1})
        bin_width           = varargin{1};
        gauss_SD            = varargin{2};
        transform           = varargin{3};    
    else
        smoothed_FR         = varargin{1};
        dim_red_FR          = varargin{2};
        empir_angle_dist_all = varargin{3};
    end
elseif nargin == 9
    smoothed_FR             = varargin{1};
    dim_red_FR              = varargin{2};
    % define the last dimension the code will look at: this assumes that we
    % don't care for all the dimensions of the neural space but only a few  
    resort_eigenv           = varargin{3};
    last_dim                = varargin{4};
end


% check that dimensions are consistent
if nargin == 7 || nargin == 9
    if ~isempty(bdf)
        if length(bdf) ~= length(smoothed_FR),error('smoothed_FR has wrong size'),end
        if length(bdf) ~= length(dim_red_FR),error('dim_red_FR has wrong size'),end
    else
        if ~isempty(smoothed_FR)
            if length(smoothed_FR) ~= length(dim_red_FR),error('dim_red_FR and smoothed_FR have different sizes'),end
        end
    end
end

% if don't have passed the option to resort the eigenvectors by
% similarity, set the flag to false
if ~exist('resort_resort_eigenv','var')
    resort_eigenv           = false;
end


% -------------------------------------------------------------------------
% Some definitions


% Nbr of BDFs (tasks)
if ~isempty(bdf)
    nbr_bdfs                = length(bdf);
else
    nbr_bdfs                = length(dim_red_FR);
end

% neural channels to be discarded
if ~isempty(bdf)
    if ~isempty(neural_chs)
        discard_neurons         = setdiff(1:length(bdf(1).units), neural_chs);
    end
end

% Nbr of neural channels
if ~isempty(neural_chs)
    nbr_neural_chs          = numel(neural_chs);
else
    if ~isempty(bdf)
        nbr_neural_chs      = length(bdf(1).units);
    else
        nbr_neural_chs      = length(dim_red_FR{1}.chs);
    end
end


% -------------------------------------------------------------------------
% Do PCA of the of the firing rates if not passed 
%
% ~-~-> optional, the function allows passing the PCA- or FA-reduced data)
%


% Smooth firing rates
if ~exist('smoothed_FR','var')
    for i = 1:nbr_bdfs
        if nargin == 5
            smoothed_FR{i}  = gaussian_smoothing2( bdf(i) ); %#ok<AGROW>
        elseif nargin == 8 && exist('bin_width','var')
            smoothed_FR{i}  = gaussian_smoothing2( bdf(i), transform, ...
                                         bin_width, gauss_SD ); %#ok<AGROW>
%             error('this option is not compatible yet with the new gaussian smoothing function');
        end
    end
end

% Do dimensionality reduction
if ~exist('dim_red_FR','var')
    for i = 1:nbr_bdfs
        dim_red_FR{i}       = dim_reduction( smoothed_FR{i}, method, ...
                                    discard_neurons ); %#ok<AGROW>
    end
end


% Check if the user has chosen to compute angles between manifolds with
% dimensionality (n) up to the dimensionality of the space (N). If that is
% the case, correct so it is up to n = N-1 as the angles when n = N are
% zero because the manifolds become the entire space
if dims_hyper_in_orig == length(dim_red_FR{1}.chs)
    dims_hyper_in_orig      = dims_hyper_in_orig - 1;
end


% -------------------------------------------------------------------------
% Load or compute distribution randomly generated manifolds, to assess the
% meaning of the angles computed below, if not passed as fcn input
%

if ~exist('empir_angle_dist_all','var')
    empir_angle_dist_all    = load(empir_angle_dist_file);
end
% check if the distribution we need is stored in empir_angle_dist_all,
% otherwise generate it
if isempty(find(empir_angle_dist_all.space_dim == nbr_neural_chs,1))
    disp('calculating random distribution 1:n-1 dimensional manifolds in n-dimensional neural space');
    [empir_angle_dist, angles_non_rand] = empirical_angle_distribution( nbr_neural_chs, ...
                                1:nbr_neural_chs, 10000 );
    empir_angle_dist.Pval   = 0.01;
else
    angles_non_rand         = empir_angle_dist_all.angle_non_orth{...
                                find(empir_angle_dist_all.space_dim == nbr_neural_chs,1)};
    empir_angle_dist         = empir_angle_dist_all.dist_angles{...
                                find(empir_angle_dist_all.space_dim == nbr_neural_chs,1)};
end


% -------------------------------------------------------------------------
% Compute angle between manifolds
%
% ~-~-> this will be done twice, using each manifold (task) as 'reference'
% to compare the angle between each of the other manifolds [see below]
%


% Find the closest eigenvector in task i+p (p>1) to each eigenvector in
% task i. Do until n = N - 1 
[angles_fc, dim_min_angle]  = find_closest_neural_hyperplane_all( dim_red_FR, ...
                                1:dims_hyper_in_orig, '' );


% Choose whether to resort the eigenvectors so they are not ranked by the
% eigenvalues of the first task but rather by similarity. Don't do by
% default
% ~-~-> likely to be deleted in newer versions of the code
if resort_eigenv
    [~, dim_min_angle]      = resort_eigenv_similarity( angles_fc, dim_min_angle, last_dim );
end

% compare hyperplanes for increasing hyperplane dimensionality, re-ordering
% the eigenvectors so as to minimize the angle between pairs according to
% dim_min_angle
angles                      = comp_neural_spaces_fcn_dim( bdf, neural_chs, dim_min_angle, ...
                                    labels, method, angles_non_rand, smoothed_FR, dim_red_FR );

                                
% Do again finding the closest eigenvector in task i-p (p>1) to each
% eigenvector in task i
[~, dim_min_angle_rev]      = find_closest_neural_hyperplane_all( dim_red_FR, ...
                                1:dims_hyper_in_orig, '', true );

angles_rev                  = comp_neural_spaces_fcn_dim( bdf, neural_chs, dim_min_angle_rev, ...
                                    labels, method, angles_non_rand, smoothed_FR, ...
                                    dim_red_FR, false, true );


                                
% -------------------------------------------------------------------------
% Define output variables

% merge the angles for the dimensionality search in both eigenvector search
% directions

angles.data                 = angles.data + angles_rev.data;

for i = 2:nbr_bdfs
    for ii = 1:i-1
        angles.labels{i,ii} = angles_rev.labels{i,ii};
    end
end


% create a struct with the empirical angle distribution and the
% orthogonality threshold
ead.data                    = empir_angle_dist;
ead.angle_non_rand          = angles_non_rand; 
ead.Pval                    = empir_angle_dist_all.Pval;
clear empir_angle_dist;
empir_angle_dist            = ead;


% look for the direction of eigenvector search (i.e., using what task
% in each pair of tasks as reference) gives the smallest manifold angle
% -- the first search criterion is find what of the two directions goes
% above the randomness threshold in a highest dimension

comb_bdfs                   = nchoosek(1:nbr_bdfs,2);
comb_bdfs_rev               = fliplr(comb_bdfs);
nbr_comb_bdfs               = size(comb_bdfs,1);

angles.min_angle            = zeros(nbr_comb_bdfs,dims_hyper_in_orig);
angles.pair_min_angle       = cell(nbr_comb_bdfs,1);

diff_w_rand_angle           = zeros(dims_hyper_in_orig,nbr_comb_bdfs);
diff_w_rand_angle_rev       = zeros(dims_hyper_in_orig,nbr_comb_bdfs);

% it happens that the angle btw the manifolds computed using both manifolds
% as reference to order the eigenvectors don't go above the randomness
% threshold
last_dim_cumsum             = 20;


for i = 1:nbr_comb_bdfs
    % compute diff w. random angle (values < 0 indicate the angle is above
    % "the randomness threshold")
    diff_w_rand_angle(:,i)  = deg2rad(angles_non_rand(1:dims_hyper_in_orig))' ...
                                 - squeeze(angles.data(comb_bdfs(i,1),comb_bdfs(i,2),...
                                1:dims_hyper_in_orig));
    diff_w_rand_angle_rev(:,i) = deg2rad(angles_non_rand(1:dims_hyper_in_orig))' ...
                                 - squeeze(angles.data(comb_bdfs_rev(i,1),comb_bdfs_rev(i,2),...
                                1:dims_hyper_in_orig));
    if find(diff_w_rand_angle(:,i)<0,1) > find(diff_w_rand_angle_rev(:,i)<0,1)
        angles.min_angle(i,:) = angles.data(comb_bdfs(i,1),comb_bdfs(i,2),...
                                1:dims_hyper_in_orig);
        angles.pair_min_angle{i} = angles.labels{comb_bdfs(i,1),comb_bdfs(i,2)};
    elseif find(diff_w_rand_angle(:,i)<0,1) < find(diff_w_rand_angle_rev(:,i)<0,1)
        angles.min_angle(i,:) = angles.data(comb_bdfs_rev(i,1),comb_bdfs_rev(i,2),...
                                1:dims_hyper_in_orig);
        angles.pair_min_angle{i} = angles.labels{comb_bdfs_rev(i,1),comb_bdfs_rev(i,2)};
    else % if neither angles or angles_rev go above rand_thr
        if squeeze(sum(angles.data(comb_bdfs(i,1),comb_bdfs(i,2),1:dims_hyper_in_orig))) > ...
                squeeze(sum(angles.data(comb_bdfs_rev(i,1),comb_bdfs_rev(i,2),1:dims_hyper_in_orig)))
            angles.min_angle(i,:) = angles.data(comb_bdfs_rev(i,1),comb_bdfs_rev(i,2),...
                                1:dims_hyper_in_orig);
            angles.pair_min_angle{i} = angles.labels{comb_bdfs_rev(i,1),comb_bdfs_rev(i,2)};
        else
            angles.min_angle(i,:) = angles.data(comb_bdfs(i,1),comb_bdfs(i,2),...
                                1:dims_hyper_in_orig);
            angles.pair_min_angle{i} = angles.labels{comb_bdfs(i,1),comb_bdfs(i,2)};
        end
    end
end


% -------------------------------------------------------------------------
% Plot the angles betweeen the manifolds in the way that they are smallest

if nbr_comb_bdfs > 1 % don't use parula for n = 1 becuse it's yellow
    cols_plot               = parula(nbr_comb_bdfs);
else
    cols_plot               = [0 0 0];
end

legends_plot                = cell(nbr_comb_bdfs,1);
for i = 1:nbr_comb_bdfs
    legends_plot{i}         = [angles.pair_min_angle{i}{1} ' vs. ' ...
                                angles.pair_min_angle{i}{2}];
end
legends_plot{length(legends_plot)+1} = ['rand. th. (P<' num2str(empir_angle_dist.Pval) ')'];

figure,hold on
for i = 1:nbr_comb_bdfs
    plot(1:dims_hyper_in_orig,rad2deg(angles.min_angle(i,:)),'color',cols_plot(i,:),...
        'linewidth',3)
end
plot(1:dims_hyper_in_orig,angles_non_rand(1:dims_hyper_in_orig),...
    'color',[.5 .5 .5],'linestyle','-.','linewidth',3)
legend(legends_plot,'Location','SouthEast','FontSize',16);legend boxoff
ylim([0 90]),xlim([0 dims_hyper_in_orig])
set(gca,'TickDir','out'), set(gca,'FontSize',16)
xlabel('manifold dimensionality'),ylabel('angle (deg)')