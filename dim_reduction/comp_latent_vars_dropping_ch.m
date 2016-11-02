%
% Compare latent variables, obtained with PCA, after dropping a certain
% percentage of channels.
%
%   CCs = comp_latent_vars_dropping_ch( bin_FRs, perc_drop, varargin )
%
%
% Inputs (opt)      : [default]
%   bin_FRs         : binned firing rates (matrix with size samples x chs)
%   mani_dim        : manifold dimensionality (= nbr. latent variables)
%   perc_drop       : percentage of channels to drop (0-1). It can be an
%                       array, in which case the function will return a 3D
%                       matrix
%   (nbr_reps)      : [100] number of combinations
%   (neural_chs)    : [all] neural channels to use
%   (normal)        : ['sqrt'] normalization of the spike firings
%
% Outputs:
%   CCs             : canonical correlations (with size: iteration x
%                       dimension x perc chs)
%
%

function CCs = comp_latent_vars_dropping_ch( bin_FRs, mani_dim, perc_drop, varargin )


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

if nargin == 6
    normal          = varargin{3};
else 
    normal          = 'sqrt';
end


% -------------------------------------------------------------------------

% get indexes all channels
all_chs             = 1:size(bin_FRs,2);
% indexes of channels to discard
chs_discard         = setdiff(all_chs,neural_chs);

% how many channels to keep
nbr_neural_chs      = numel(neural_chs);
nbr_chs_keep        = round(nbr_neural_chs*(1-perc_drop));

% square root transform the firing rates, if specified
switch normal
    case 'sqrt'
        bin_FRs     = sqrt(bin_FRs);
    otherwise
end

% % temporary: keep more latent variables when dropping channels
% mani_dim_drop       = 20;
mani_dim_drop       = mani_dim;

% -------------------------------------------------------------------------
% Do

% preallocate matrices
chs                 = zeros(nbr_reps,max(nbr_chs_keep));
CCs                 = zeros(nbr_reps,mani_dim,length(perc_drop));

% add a row of zeros to bin_FRs, for compatibility with dim_reduction.m
bin_FRs             = [zeros(size(bin_FRs,1),1), bin_FRs];

% get leading 'mani_dim' latent variables for all the channels
dim_red_all         = dim_reduction( bin_FRs, 'pca', chs_discard );

for p = 1:length(perc_drop)
    for i = 1:nbr_reps

        % randomly choose a subset with 'perc_chs' % of the channels
        chs_this            = datasample( 1:nbr_neural_chs, nbr_chs_keep(p), ...
                                'Replace', false );

        % do PCA on the channels specified in 'chs'
        chs_discard_this    = setdiff(all_chs,chs_this);
        dim_red_this        = dim_reduction( bin_FRs, 'pca', chs_discard_this );

        % do CC between the leading 'mani_dim' latent variables computed from
        % all channels and the leading 'mani_dim' latent variables computed
        % from the subset of channels
        [~,~,CCs(i,:,p),~,~,stats_CC{i,p}] = canoncorr( dim_red_all.scores(:,1:mani_dim),...
                            dim_red_this.scores(:,1:mani_dim_drop) );
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

figure, hold on
for p = 1:length(perc_drop)
    plot(mean(CCs(:,:,p),1),'linewidth',4,'color',colors(p,:))
end
legend(leg,'Location','Southwest'), legend boxoff
for p = 1:length(perc_drop)
    plot(mean(CCs(:,:,p),1)+std(CCs(:,:,p),0,1),'linestyle','-.','linewidth',2,'color',colors(p,:))
    plot(mean(CCs(:,:,p),1)-std(CCs(:,:,p),0,1),'linestyle','-.','linewidth',2,'color',colors(p,:))
end
set(gca,'TickDir','out'), set(gca,'FontSize',14)
xlabel('latent variable'),ylabel('canonical correlation')
xlim([0 mani_dim+1]),ylim([0 1])
box off

drawnow