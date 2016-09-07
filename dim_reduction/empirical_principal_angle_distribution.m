%
% Create a random distribution of angles between n-dimensional hyperplanes
% in an m-dimensional hyperspace (n<=m)
%
%   function [dist_angles angle_non_orth] = empirical_principal_angle_distribution( space_dim, plane_dim, samples )
%
% Inputs (opt)      : [default]
%   space_dim       : dimensionality of the space. Can be a scalar or
%                       a vector; if it is a scalar, the function will
%                       create a distribution of angles between hyperplanes
%                       in spaces with each of these dimensionalities, if
%                       it is a vector, it will create one distribution of
%                       'samples' samples for each of the combination of in
%                       plane_dim x space_dime
%   plane_dim       : dimensionality of the hyperplanes. Can be a scalar or
%                       a vector; if it is a scalar, the function will
%                       create a distribution of angles between hyperplanes
%                       with that dimensionality, if it is a vector, it
%                       will create one distribution of 'samples' samples
%                       for each of the dimensionalities in plane_dim
%   samples         : number of angles in the distribution
%   (P_orth)        : [0.001] area under the PDF above which we'll consider
%                       the hyperplanes not to be different from orthogonal
%   (plot_yn)       : [false] plot the PDF (bool)
%
% Outputs:
%   dist_princ_angles : vector or matrix with the distribution of principal
%                       angles 
%   angle_non_orth  : maximum angle at which P < P_orth
%
%

function [dist_princ_angles, angle_non_orth] = empirical_principal_angle_distribution( ...
                                            space_dim, plane_dim, samples, varargin )


% input parameters
if nargin >= 4
    if ~isempty(varargin{1}), P_orth = varargin{1}; else
    P_orth          = 0.001; end
else
    P_orth          = 0.001;
end
if nargin == 5
    plot_yn         = varargin{2};
else
    plot_yn         = false;
end

% preallocate to store results
% dist_princ_angles has size number of samples x nbr of principal angles x
% dimensionality of the hyperplanes x dimensionality of the hyperspaces
dist_princ_angles   = zeros(samples, max(plane_dim), length(plane_dim), length(space_dim));


% 1. create random planes and compute their dist_angles. Note that the
% planes don't need to be defined by unitary vectors (i.e. the columns of
% matrices A and B don't need to be unitary vectors), as 'subspace'
% orthogonalizes the matrices 
for s = 1:length(space_dim)
    disp(['generating hyperplane distribution for hyperspace with dimensionality: ' num2str(space_dim(s))]);
    for p = 1:length(plane_dim)
        disp([' and hyperplanes with dimensionality: ' num2str(plane_dim(p))]);
        disp('...');
        for i = 1:samples
            % create planes
            A       = randn(space_dim(s),plane_dim(p));
            B       = randn(space_dim(s),plane_dim(p));
            dist_princ_angles(i,1:plane_dim(p),p,s) = principal_angles(A,B);
        end
    end
    disp(' ');
end


% -------------------------------------------------------------------------
% turn the angles into a distribution of angles
hist_x              = 0:pi/2/900:pi/2;

% hist_dist_princ_angles has size number of hist bins x nbr of principal
% angles (=dimensionality of the hyperplane) x dimensionality of the
% hyperplanes x dimensionality of the hyperspaces 
hist_dist_princ_angles = nan(length(hist_x)-1,max(plane_dim),length(plane_dim),length(space_dim));
for s = 1:length(space_dim)
    for p = 1:length(plane_dim)
        % get number of dimensions of the hyperplanes for this s and p --note
        % that they can be different if specified
        nbr_dims_this = numel(find(dist_princ_angles(1,:,p,s)));
        for d = 1:nbr_dims_this
            hist_dist_princ_angles(:,d,p,s) = histcounts(dist_princ_angles(:,d,p,s), hist_x)/samples;
        end
    end
end
    
% find the maximum angle at which P < 0.01
angle_non_orth      = zeros(max(plane_dim),length(plane_dim),length(space_dim));
for s = 1:length(space_dim)
    for p = 1:length(plane_dim)
        % get number of dimensions of the hyperplanes for this s and p --note
        % that they can be different if specified
        nbr_dims_this = numel(find(dist_princ_angles(1,:,p,s)));
        for d = 1:nbr_dims_this 
            angle_non_orth(d,p,s) = hist_x(find(cumsum(hist_dist_princ_angles(:,d,p,s))>P_orth, 1));
        end
    end
end
% turn into degrees
% note: the histogram starts in 1 deg and hist_x in 0, which means we'd
% have to add 1 to max_angle_non_orth, but we compensate for that by
% looking for the first bin at which P>P_orth
angle_non_orth      = rad2deg(angle_non_orth); 


% plot distributions ---very messy, needs to be fixed
if plot_yn
    
    % plot the random principal angles, for each hyperspace dimensionality 
    for p = 1:length(plane_dim)
        aux_colors  = parula(length(space_dim));
        figure, hold on
        for s = 1:length(space_dim)
            plot(squeeze(angle_non_orth(:,p,s)),'linewidth',2,'color',aux_colors(s,:))
        end
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        xlabel('dimension'),ylabel(['min angle P < ' num2str(P_orth) ' btw ' ...
            num2str(plane_dim(p)) ' hyperplanes'])
        xlim([0 plane_dim(p) + 1]),ylim([0 90])
        aux_legend = cell(1,length(space_dim));
        for s = 1:length(space_dim)
            aux_legend{s} = ['n = ' num2str(space_dim(s))];
        end
        legend(aux_legend,'Location','Northwest'), legend boxoff
    end
    
    % plot the distributions, for each principal angle, for one of the
    % plane and space dimensionalities
    for s = 1:length(space_dim)
        
        this_plane_dim = 1;
        this_space_dim = s;
        nbr_dims_this = numel(find(dist_princ_angles(1,:,this_plane_dim,this_space_dim)));
        aux_colors_hist = parula(nbr_dims_this);
        figure, hold on
        for d = 1:nbr_dims_this
            plot(rad2deg(hist_x(1:(length(hist_x)-1))),hist_dist_princ_angles(:,d,this_plane_dim,this_space_dim),...
                'color',aux_colors_hist(d,:))
        end
        xlim([0 90]), xlabel('angle (deg)'),
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        ylabel(['normalized counts - space dim: ' num2str(space_dim(s))])
        aux_legend = cell(1,nbr_dims_this);
        for d = 1:nbr_dims_this
            aux_legend{d} = ['p.a. = ' num2str(d)];
        end
        legend(aux_legend,'Location','Northwest'), legend boxoff
    end
    
    
    % more plotting stuff to be implemented....
end



% % -------------------------------------------------------------------------
% % code to plot the principal angles for a few space dimensionalities
% 
% % space number
% indx_plane_dim  = 1;
% indx_space_dims = 1:size(dist_princ_angles,4);
% 
% last_plane_dim = 20;
% 
% colors  = [ .6 .6 .6    % grey
%             1 .6 0      % orange
%             .7 .9 .9    % light blue
%             .6 1 .6     % light green
%             .9 .9 .5];  % yellowish
%         
% 
% figure,hold on
% for s = 1:length(indx_space_dims)
%     for n = 1:1000,
%         plot(rad2deg(dist_princ_angles(n,1:last_plane_dim,s)),'color',colors(s,:)) 
%     end
% end
% set(gca,'Tickdir','out'),set(gca,'FontSize',14)
% ylim([0 90])