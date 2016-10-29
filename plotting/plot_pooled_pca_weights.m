%
% Function that plots PCA weights for all tasks in a 'datasets' struct
%
%   plot_pooled_pca_weights( datasets )    
%   
%
% Inputs (opt)          : [default]
%   datasets            : cell array of 'datasets' structs. Can also be a
%                           cell of dim_red structs
%   (dims)              : [1:20] manifold dimensions
%
%


function plot_pooled_pca_weights( datasets, varargin )


% read inputs
if nargin == 2
    dims        = varargin{1};
else
    dims        = 1:20;
end


% create a cell array h with the PC
for i = 1:length(datasets)
    % check if it is a datasets cell
    if iscell(datasets)
        h{i}.w  = cellfun( @(x) x.w(:,dims), datasets{i}.dim_red_FR, 'UniformOutput', false );
    % or a dim_red_FR cell
    else
        h{1}.w  = cellfun( @(x) x.w(:,dims), datasets.dim_red_FR, 'UniformOutput', false );
    end
end


% color map for plotting
nbr_tasks       = sum( cellfun( @(x) length(x.w), h ) );
colors          = parula( nbr_tasks );


% do
figure,hold on, c = 1;
for i = 1:length(h)
    for j = 1:length(h{i}.w)
        % get histogram
        h_x     = -1.01:0.01:1; 
        hst = histcounts(h{i}.w{j},h_x);
        plot(h_x(2:end),hst/sum(hst),'color',colors(c,:),'linewidth',1.5)
        c       = c + 1;
    end
end
set(gca,'FontSize',14,'TickDir','out')
xlabel('PCA weight')
ylabel('norm. counts')
