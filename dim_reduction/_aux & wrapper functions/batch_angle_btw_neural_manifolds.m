%
% Batch analysis to compute the angles between all pairs of tasks that
% monkeys performed in the same session. The function will run for all the
% datasets in "all_manifold_datasets.mat"
%
%
%

function angle_results = batch_angle_btw_neural_manifolds( varargin )


% -------------------------------------------------------------------------
% some options that could be turn into fcn params
dim_manifold            = 20; % the neural manifolds are defined by the first 20 dimensions
last_eigenv             = 20; % when counting things based on eigenvectors go up to 20

empir_angle_dist_file   = ['/Users/juangallego/Documents/NeuroPlast/Data/' ...
                            '_Dimensionality reduction/_control analyses/' ...
                            'empirical angle distribution all datasets.mat'];
Pval_empir_angle_dist   = 0.001; % can be 0.01 --this needs to be improved

% -------------------------------------------------------------------------
% read input params

if nargin
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end
clear varargin;


% -------------------------------------------------------------------------
% load data

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end

empir_angle_dist_all    = load(empir_angle_dist_file);

% choose the random distributions based on the threshold 
% -- ToDo: improve this code
switch Pval_empir_angle_dist
    case 0.01
        % do not need to do anything because this is the by-default Pval
    case 0.001
        % replace the angle_non_orth, which is for P = 0.01, for the
        % distribution for P = 0.001
        empir_angle_dist_all.angle_non_orth_01 = empir_angle_dist_all.angle_non_orth;
        empir_angle_dist_all.angle_non_orth = empir_angle_dist_all.angle_non_orth_001;
end
% add Pval to empir_angle_dist_all for plotting
empir_angle_dist_all.Pval = Pval_empir_angle_dist;


% -------------------------------------------------------------------------
% get info about monkeys, tasks, and who did what

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% compute angles between manifoldsclc

for i = 1:meta_info.nbr_monkeys
    
    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
        
        % what dataset are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);
        
        % -----------------------------------------------------------------
        % 1) angles between pairs of eigenvectors --show angles and
        % corresponding eigenvectors with each of the tasks as reference
        
        % get non-orthogonality angle
        angle_orth      = empir_angle_dist_all.angle_non_orth_001{ ...
            find( empir_angle_dist_all.space_dim == ...
            length(datasets{dtst}.dim_red_FR{1}.chs) )}(1);
        
        % compute angles btw eigenvectors
        angles_eigenv   = comp_indiv_dims_neural_manifolds( ...
            datasets{dtst}.dim_red_FR, 1:dim_manifold, datasets{dtst}.labels, ...
            angle_orth );
        
        % -----------------------------------------------------------------
        % 2) angles between manifolds
        [angles, ~, ~, empir_angle_dist] = comp_neural_manifolds( ...
            [], datasets{dtst}.neural_chs, dim_manifold, datasets{dtst}.labels, '', ...
            [], datasets{dtst}.dim_red_FR, empir_angle_dist_all );

        % -----------------------------------------------------------------
        % 3) store all results
        
        % for the manifold comparison
        data{dtst}.angles_manifolds = angles;
        data{dtst}.angle_non_orth_manifolds = empir_angle_dist.angle_non_rand(1:dim_manifold);
        
        % for the eigenvector comparison 
        data{dtst}.angles_eigenv = angles_eigenv;
        data{dtst}.angle_non_orth_eigenv = angle_orth;
         
    end
end


% -------------------------------------------------------------------------
% Do some summary analyses

% -------------------------------------------
% 0. some preliminary definitions:

% how many pairs of combinations ( monkey x session x pairs of tasks ) do we have?
nbr_manifold_pairs      = length(meta_info.task_pairs.monkey);

% ---------------------
% For the manifold analysis
% var to store the dimensionality of the highest dimensional manifold below
% the "randomness threshold" for each pair of manifolds --note that it will
% be a NaN if it doesn't go above the threshold for the 1:dim_manifold
% first dimensions
last_dim_below_th       = nan(1,nbr_manifold_pairs);
manifold_pair           = cell(1,nbr_manifold_pairs);

% ---------------------
% For the eigenvector analysis
nbr_non_orth_eigenv     = zeros(1,nbr_manifold_pairs);
task_pair_eigenv_search = cell(1,nbr_manifold_pairs);

% -------------------------------------------
% 1. find: 
%     a) the dimension at which the angle goes above the random threshold
%     for all pairs of tasks, and 
%     b) the number of eigenvectors that are non orthogonal (1:last_eigenv)

ctr                     = 1;
for d = 1:length(data)
    for p = 1:size(data{d}.angles_manifolds.min_angle,1)
        
        % a) find what manifold dimension goes above the randomness
        % threshold
        dim_above_th    = find( rad2deg(data{d}.angles_manifolds.min_angle(p,:)) > data{d}.angle_non_orth_manifolds, 1);
        if ~isempty(dim_above_th)
            last_dim_below_th(ctr) = dim_above_th - 1;
        end
        manifold_pair{ctr} = data{d}.angles_manifolds.pair_min_angle{p};
        
        % b) the number of eigenvectors that are non orthogonal

        % -------------
        % TODO: THIS HAS TO BE FINISHED
        
        % task_pair_eigenv_search{ctr} = data{d}.angle_eigenv
        
        % UP THERE
        % -------------
        
        % update ctr
        ctr             = ctr + 1;
    end
end


% -------------
% AFTER THIS; EVERYTHING FOR THE MANIFOLDS HAS TO BE REPLICATED FOR THE
% EIGENVECTORS...
% -------------


% do a histogram to summarize this
[counts_dim_below_th, edges_dim_below_th] = histcounts(last_dim_below_th,1:dim_manifold+1);
edges_dim_below_th      = edges_dim_below_th(1:dim_manifold);

counts_non_orth_pairs   = sum(isnan(last_dim_below_th));

% store for later
summary_data.hist_all.counts_dim_below_th = counts_dim_below_th;
summary_data.hist_all.edges_dim_below_th = edges_dim_below_th;
summary_data.hist_all.counts_non_orth_pairs = counts_non_orth_pairs;

% -------------------------------------------
% 2. the same, keeping track of the monkey and pair of tasks

% 'last_dim_below_th' is ordered the same way as
% meta_info.task_pairs.task_pair. We can use that for assigning the monkey
% and task

summary_data.manifold_pair_min_angle = manifold_pair;
summary_data.last_dim_below_th = last_dim_below_th;
summary_data.per_pair   = struct('pair',[meta_info.task_pairs.unique_pairs],...
                            'last_dim_below_th',[],'pair_nbr',[]);

for i = 1:numel(summary_data.manifold_pair_min_angle)
    for ii = 1:numel(meta_info.task_pairs.unique_pairs)
        pair_cmp        = strcmpi(meta_info.task_pairs.unique_pairs{ii},...
                            sort(summary_data.manifold_pair_min_angle{i}));
    	if sum(pair_cmp) == 2
            summary_data.per_pair(ii).last_dim_below_th = [...
                summary_data.per_pair(ii).last_dim_below_th, ...
                summary_data.last_dim_below_th(i)];
            summary_data.per_pair(ii).pair_nbr = [...
                summary_data.per_pair(ii).pair_nbr, i];
            break;
        end
    end
end

% now do histograms for each pair of tasks
for p = 1:numel(meta_info.task_pairs.unique_pairs)
    [counts_this_pair, edges_this_pair] = histcounts(...
        summary_data.per_pair(p).last_dim_below_th,1:dim_manifold+1);
    edges_this_pair         = edges_this_pair(1:dim_manifold);
    % store
    summary_data.per_pair(p).hist.counts_below_th = counts_this_pair;
    summary_data.per_pair(p).hist.edges = edges_this_pair;
    summary_data.per_pair(p).hist.counts_non_orth = sum(isnan(...
        summary_data.per_pair(p).last_dim_below_th));
end



% add variables to results

angle_results.data     	= data;
angle_results.summary_data = summary_data;


% -------------------------------------------------------------------------
% PLOTS

% 1. histogram with the dimensionality of the highest dimensional manifold
% with which angle is below the threshold
figure, hold on
bar(edges_dim_below_th,counts_dim_below_th,'Facecolor','k','Edgecolor','k')
set(gca,'TickDir','out'), set(gca,'FontSize',16)
bar(dim_manifold+3,counts_non_orth_pairs,'Facecolor',[.6 .6 .6],'Edgecolor',[.6 .6 .6])
ylim([0 max([counts_non_orth_pairs, max(counts_dim_below_th)])+1]), xlim([0 dim_manifold+5])
ylabel('counts'),xlabel(['highest dimensionality non-orthogonal manifolds (P<' num2str(Pval_empir_angle_dist) ')'])
set(gca,'XTick',[0:5:dim_manifold dim_manifold+3])
set(gca,'XTickLabel',[num2cell([0:5:20]) '>20'])


% 2. Same histogram but stacking according to the task

% matrix for the stack plot
matrix_stack = zeros(dim_manifold,numel(meta_info.task_pairs.unique_pairs));
for p = 1:numel(meta_info.task_pairs.unique_pairs)
    matrix_stack(:,p)       = summary_data.per_pair(p).hist.counts_below_th;
end

matrix_stack_non_orth = zeros(1,numel(meta_info.task_pairs.unique_pairs));
for p = 1:numel(meta_info.task_pairs.unique_pairs)
    matrix_stack_non_orth(:,p) = summary_data.per_pair(p).hist.counts_non_orth;
end
% dirty fix to make the stack command work
matrix_stack_non_orth       = [matrix_stack_non_orth; zeros(1,length(matrix_stack_non_orth))];

% legend
legend_plot                 = cell(numel(meta_info.task_pairs.unique_pairs),1);
for i = 1:numel(meta_info.task_pairs.unique_pairs)
    legend_plot{i}          = [summary_data.per_pair(i).pair{1} ' vs. ' ...
                                summary_data.per_pair(i).pair{2}];
end

figure,hold on
bar(matrix_stack,'stack','EdgeColor',[.6 .6 .6])
bar([dim_manifold+3,dim_manifold+4],matrix_stack_non_orth,'stack','Edgecolor',[.6 .6 .6])
ylim([0 max([counts_non_orth_pairs, max(counts_dim_below_th)])+1]), xlim([0 dim_manifold+5])
ylabel('counts'),xlabel(['highest dimensionality non-orthogonal manifolds (P<' num2str(Pval_empir_angle_dist) ')'])
set(gca,'XTick',[0:5:dim_manifold dim_manifold+3])
set(gca,'XTickLabel',[num2cell([0:5:20]) '>20'])
set(gca,'TickDir','out'), set(gca,'FontSize',16)
legend(legend_plot,'Location','NorthWest'),legend boxoff


% Histogram with the number of non-orthogonal pairs of eigenvectors per
% type of task
