%
% Compare neural manifolds, obtained with PCA, taking pairs subsets with a
% certain % of all the channels
%
%   PAs = comp_manifolds_dropping_ch( bin_FRs, perc_drop, varargin )
%
%
% Inputs (opt)      : [default]
%   bin_FRs         : binned firing rates (matrix with size samples x chs)
%   mani_dim        : manifold dimensionality (= nbr. latent variables)
%   perc_drop       : percentage of channels to drop (0-1). It can be an
%                       array, in which case the function will return a 3D
%                       matrix
%   (nbr_reps)      : [100] number of combinations
%   (neural_chs)    : [all] neural channels to use. If empty it will use
%                       all
%   (normal)        : ['sqrt'] normalization of the spike firings
%   (ead)            : empirical angle distribution structure
%
% Outputs:
%   PAs             : principal angles (with size: iteration x
%                       dimension x perc chs)
%
%


function PAs = comp_manifolds_dropping_ch( bin_FRs, mani_dim, perc_drop, varargin )


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
    neural_chs      = varargin{2};
    % if empty, choose all channels
    if isempty(neural_chs)
        neural_chs      = 1:size(bin_FRs,2);
    end
else
    neural_chs      = 1:size(bin_FRs,2);
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


% -------------------------------------------------------------------------

% get indexes all channels
all_chs             = 1:size(bin_FRs,2);

% how many channels to keep
nbr_neural_chs      = numel(neural_chs);
nbr_chs_keep        = round(nbr_neural_chs*(1-perc_drop));

% square root transform the firing rates, if specified
switch normal
    case 'sqrt'
        bin_FRs     = sqrt(bin_FRs);
    otherwise
end

% retrieve random angles
if ~isempty(find(ead.space_dim==nbr_neural_chs))
    if ead.plane_dim == mani_dim
        ang_rand    = ead.angle_non_orth(:,1,find(ead.space_dim==nbr_neural_chs));
    else
        warning('do not have empirical angle distribution for this manifold dimensionality');
    end
else
    warning('do not have empirical angle distribution for this space dimensionality');
end

% -------------------------------------------------------------------------
% Do

% preallocate matrices
chs                 = zeros(nbr_reps,max(nbr_chs_keep));
PAs                 = zeros(nbr_reps,mani_dim,length(perc_drop));

% add a row of zeros to bin_FRs, for compatibility with dim_reduction.m
bin_FRs             = [zeros(size(bin_FRs,1),1), bin_FRs];


for p = 1:length(perc_drop)
    for i = 1:nbr_reps

        % randomly choose a subset with 'perc_chs' % of the channels
        chs_1               = datasample( 1:nbr_neural_chs, nbr_chs_keep(p), ...
                                'Replace', false );
        chs_2               = datasample( 1:nbr_neural_chs, nbr_chs_keep(p), ...
                                'Replace', false );

        % do PCA on the channels specified in 'chs'
        chs_discard_1       = setdiff(all_chs,chs_1);
        dim_red_1           = dim_reduction( bin_FRs, 'pca', chs_discard_1 );

        chs_discard_2       = setdiff(all_chs,chs_2);
        dim_red_2           = dim_reduction( bin_FRs, 'pca', chs_discard_2 );
        
        % do CC between the leading 'mani_dim' latent variables computed from
        % all channels and the leading 'mani_dim' latent variables computed
        % from the subset of channels
        PAs(i,:,p)          = rad2deg(principal_angles( dim_red_1.w(:,1:mani_dim),...
                                dim_red_2.w(:,1:mani_dim) ));
    end
end


% -------------------------------------------------------------------------
% summary plot
if length(perc_drop) > 1
    colors              = parula(length(perc_drop));
else
    colors              = [0 0 0];
end
for i = 1:length(perc_drop)
    leg{i}              = [num2str(perc_drop(i)) '% ch. drop'];
end
if exist('ang_rand','var')
    leg{length(leg)+1}  = ['P < ' num2str(ead.P_orth)];
end

figure, hold on
for p = 1:length(perc_drop)
    plot(mean(PAs(:,:,p),1),'linewidth',4,'color',colors(p,:))
end
if exist('ang_rand','var')
    plot(ang_rand,'linewidth',4,'color',[.6 .6 .6],'linestyle','-.')
end
legend(leg,'Location','Northwest'), legend boxoff
for p = 1:length(perc_drop)
    plot(mean(PAs(:,:,p),1)+std(PAs(:,:,p),0,1),'linestyle','-.','linewidth',2,'color',colors(p,:))
    plot(mean(PAs(:,:,p),1)-std(PAs(:,:,p),0,1),'linestyle','-.','linewidth',2,'color',colors(p,:))
end
set(gca,'TickDir','out'), set(gca,'FontSize',14)
xlabel('latent variable'),ylabel('principal angles')
xlim([0 mani_dim+1]),ylim([0 90])
box off

drawnow
