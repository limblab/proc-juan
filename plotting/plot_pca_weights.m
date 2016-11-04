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

% make C axis a little bit "wider" than what we need and use the same for
% all panels
c_axis_max          = ceil( max( cellfun( @(x) max(max(abs(x.w))), dim_red_FR ) )*10 )/ 10;
x_ticks             = 10:10:(floor(nbr_chs/10)*10);

% plot
figure
for i = 1:nbr_bdfs
    subplot(nbr_rows,nbr_cols,i)
    imagesc(1:nbr_chs,1:nbr_chs,abs(dim_red_FR{i}.w)),colorbar
    caxis([0 c_axis_max])
    set(gca,'TickDir','out','FontSize',14)
    set(gca,'XTick',x_ticks)
    title(labels(i)),xlabel('principal component'),ylabel('neural channel')
end