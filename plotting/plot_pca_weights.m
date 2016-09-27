%
% Plot PCA weights for all the tasks in a session
%

function plot_pca_weights( dim_red_FR, varargin )


% read input params
if nargin == 2
    labels          = varargin{1};
end
    

% get nbr BDFs
nbr_bdfs            = length(dim_red_FR);

% if task labels haven't been passed, create a generic one
if ~exist('labels','var')
    for i = 1:nbr_bdfs
        labels{i}   = ['task ' num2str(i)];
    end
end


% get some info for the plot
nbr_chs             = length(dim_red_FR{1}.eigen);

nbr_rows            = floor(sqrt(nbr_bdfs));
nbr_cols            = ceil(nbr_bdfs/nbr_rows);

% plot
figure
for i = 1:nbr_bdfs
    subplot(nbr_rows,nbr_cols,i)
    imagesc(1:nbr_chs,1:nbr_chs,abs(dim_red_FR{i}.w)),colorbar
    set(gca,'TickDir','out','FontSize',14)
    title(labels(i)),xlabel('principal component'),ylabel('neural channel')
end