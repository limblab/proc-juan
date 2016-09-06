%
% Batch analysis to compute the angles between all pairs of tasks that
% monkeys performed in the same session. The function will run for all the
% datasets in "all_manifold_datasets.mat"
%
%
% Note: a lot of the code here needs to be improved
%
%

function angle_results = batch_angle_btw_manifold_dims( varargin )


% -------------------------------------------------------------------------
% some options that could be turn into fcn params

% the dimensionality of the neural manifolds, and how many eigenvectors
% will be compared in the eigenvector comparisons
dim_manifold            = 20; 

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
% -------------------------------------------------------------------------
% compute angles between pairs of eigenvectors

data                    = cell(1,length(meta_info.tasks_per_session));

for i = 1:meta_info.nbr_monkeys
    
    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
        
        % what dataset are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);
        
        % -----------------------------------------------------------------
        % 1) angles between pairs of task-manifold dimensions 
         
        % get non-orthogonality angle for 1D vectors in the n-dimensional space
        angle_orth      = empir_angle_dist_all.angle_non_orth_001{ ...
            find( empir_angle_dist_all.space_dim == ...
            length(datasets{dtst}.dim_red_FR{1}.chs) )}(1);

        % compute principal angles between all pairs of tasks in this
        % session
        princ_angles    = principal_angles_all_manifolds( datasets{dtst}.dim_red_FR, ...
                            1:dim_manifold, datasets{dtst}.labels, angle_orth );
        
%         % -----------------------------------------------------------------
%         % 2) angles between manifolds
%         [angles, ~, ~, empir_angle_dist] = comp_neural_manifolds( ...
%             [], datasets{dtst}.neural_chs, dim_manifold, datasets{dtst}.labels, '', ...
%             [], datasets{dtst}.dim_red_FR, empir_angle_dist_all );

        % -----------------------------------------------------------------
        % 3) store all results
        
        % for the eigenvector (PC) comparison
        data{dtst}.princ_angles = princ_angles;
    end
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Do some summary analyses

% -------------------------------------------
% 0. some preliminary definitions:

% how many pairs of combinations ( monkey x session x pairs of tasks ) do we have?
nbr_manifold_pairs      = length(meta_info.task_pairs.monkey);
nbr_pairs_tasks      	= length(meta_info.task_pairs.unique_pairs);

% % ---------------------
% % For the manifold analysis
% % var to store the dimensionality of the highest dimensional manifold below
% % the "randomness threshold" for each pair of manifolds --note that it will
% % be a NaN if it doesn't go above the threshold for the 1:dim_manifold
% % first dimensions
% last_dim_below_th       = nan(1,nbr_manifold_pairs);
% manifold_pair           = cell(1,nbr_manifold_pairs);



% -------------------------------------------
% 1. find: 
%     a) the dimension at which the angle goes above the random threshold
%     for all pairs of tasks, and 
%     b) the number of eigenvectors that are non orthogonal (1:dim_manifold)

nbr_eigenv_below_thr    = zeros(1,size(meta_info.task_pairs.task_pair,2));
angles_matrix           = zeros(nbr_manifold_pairs,dim_manifold);
ctr                     = 1;
for d = 1:length(data)
    for p = 1:size(data{d}.princ_angles.angles,1)
        
%         % a) find what manifold dimension goes above the randomness
%         % threshold
%         dim_above_th    = find( rad2deg(data{d}.angles_manifolds.min_angle(p,:)) > data{d}.angle_non_orth_manifolds, 1);
%         if ~isempty(dim_above_th)
%             last_dim_below_th(ctr) = dim_above_th - 1;
%         end
%         manifold_pair{ctr} = data{d}.angles_manifolds.pair_min_angle{p};
        

        % b) the number of eigenvectors that are non orthogonal

        nbr_eigenv_below_thr(ctr)   = sum(rad2deg(data{d}.princ_angles.angles(p,:)) < ...
                                        data{d}.princ_angles.angle_orth);
        angles_matrix(ctr,:)        = data{d}.princ_angles.angles(p,:);
        % update ctr
        ctr             = ctr + 1;
    end
end


% -------------------------------------------------------------------------
% structure to store the summary data 

% create 
summary_data.princ_angles   = struct('pair',[meta_info.task_pairs.unique_pairs],...
                            'nbr_eigenv_below_thr',[],'angles',zeros(1,dim_manifold));




% -------------------------------------------------------------------------
% SUMMARIZE RESULTS FOR THE PRINCIPAL ANGLES
% -------------------------------------------------------------------------


% -- Note: 'nbr_eigenvals_above_thr' is ordered the same way as
% meta_info.task_pairs.task_pair. We can use that for assigning the monkey 
% and task


% summarize, for each pair of tasks (e.g., tube vs ball, or iso vs mov) how
% many eigenvectors are below the 'randomness threshold' (it will be a
% number per session)
for i = 1:size(nbr_eigenv_below_thr,2)
    for ii = 1:nbr_pairs_tasks
        % retrieve pairs ot tasks
        pair_cmp        = strcmpi(meta_info.task_pairs.unique_pairs{ii},...
                            sort(meta_info.task_pairs.task_pair{i}));
    	if sum(pair_cmp) == 2
            % store angles in angles matrix -- if the matrix has zeroes fill
            % it with the number
            % ToDo: improve
            if numel(unique(summary_data.princ_angles(ii).angles)) == 1
                summary_data.princ_angles(ii).angles(1,:) = angles_matrix(i,:);
            else
                summary_data.princ_angles(ii).angles = [summary_data.princ_angles(ii).angles;
                    angles_matrix(i,:)];
            end
            % get nbr canonical angles smaller than the randomness
            % threshold
            summary_data.princ_angles(ii).nbr_eigenv_below_thr = [...
                summary_data.princ_angles(ii).nbr_eigenv_below_thr, ...
                nbr_eigenv_below_thr(i)];

            break;
        end
    end
end

% now do the histogram of this
for p = 1:numel(meta_info.task_pairs.unique_pairs)
    [counts_this_pair, edges_this_pair] = histcounts(...
        summary_data.princ_angles(p).nbr_eigenv_below_thr,1:dim_manifold+1);
    edges_this_pair         = edges_this_pair(1:dim_manifold);
    % store
    summary_data.princ_angles(p).hist.counts_below_th = counts_this_pair;
    summary_data.princ_angles(p).hist.edges = edges_this_pair;
    summary_data.princ_angles(p).hist.counts_non_orth = sum(isnan(...
        summary_data.princ_angles(p).nbr_eigenv_below_thr));
end






% =========================================================================
% =========================================================================
% =========================================================================
% AQUI




% % -------------------------------------------------------------------------
% % SUMMARIZE RESULTS FOR THE MANIFOLDS
% % -------------------------------------------------------------------------
% 
% % do a histogram to summarize this
% [counts_dim_below_th, edges_dim_below_th] = histcounts(last_dim_below_th,1:dim_manifold+1);
% edges_dim_below_th      = edges_dim_below_th(1:dim_manifold);
% 
% counts_non_orth_pairs   = sum(isnan(last_dim_below_th));
% 
% % store for later
% summary_data.hist_all.counts_dim_below_th = counts_dim_below_th;
% summary_data.hist_all.edges_dim_below_th = edges_dim_below_th;
% summary_data.hist_all.counts_non_orth_pairs = counts_non_orth_pairs;
% 
% % -------------------------------------------
% % 2. the same, keeping track of the monkey and pair of tasks
% 
% % 'last_dim_below_th' is ordered the same way as
% % meta_info.task_pairs.task_pair. We can use that for assigning the monkey
% % and task
% 
% summary_data.manifold_pair_min_angle = manifold_pair;
% summary_data.last_dim_below_th = last_dim_below_th;
% summary_data.per_pair   = struct('pair',[meta_info.task_pairs.unique_pairs],...
%                             'last_dim_below_th',[],'pair_nbr',[]);
% 
% for i = 1:numel(summary_data.manifold_pair_min_angle)
%     for ii = 1:numel(meta_info.task_pairs.unique_pairs)
%         pair_cmp        = strcmpi(meta_info.task_pairs.unique_pairs{ii},...
%                             sort(summary_data.manifold_pair_min_angle{i}));
%     	if sum(pair_cmp) == 2
%             summary_data.per_pair(ii).last_dim_below_th = [...
%                 summary_data.per_pair(ii).last_dim_below_th, ...
%                 summary_data.last_dim_below_th(i)];
%             summary_data.per_pair(ii).pair_nbr = [...
%                 summary_data.per_pair(ii).pair_nbr, i];
%             break;
%         end
%     end
% end
% 
% % now do histograms for each pair of tasks
% for p = 1:numel(meta_info.task_pairs.unique_pairs)
%     [counts_this_pair, edges_this_pair] = histcounts(...
%         summary_data.per_pair(p).last_dim_below_th,1:dim_manifold+1);
%     edges_this_pair         = edges_this_pair(1:dim_manifold);
%     % store
%     summary_data.per_pair(p).hist.counts_below_th = counts_this_pair;
%     summary_data.per_pair(p).hist.edges = edges_this_pair;
% %     summary_data.per_pair(p).hist.counts_non_orth = sum(isnan(...
% %         summary_data.per_pair(p).last_dim_below_th));
% end




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% add variables to results

angle_results.data     	= data;
angle_results.summary_data = summary_data;




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOTS


% Canonical angles for all the tasks, marking the last dimension below
% chance
colors                      = parula(nbr_pairs_tasks);

figure,hold on
for i = 2:nbr_pairs_tasks
    these_angles            = rad2deg(summary_data.princ_angles(i).angles');
    these_last_eig_below_thr = summary_data.princ_angles(i).nbr_eigenv_below_thr;
    
    % plot the traces
	plot(these_angles,'color',colors(i,:),'linewidth',2)
    % and add markers showing the last dim below the randomness threshold
    for ii = 1:size(these_angles,2)
        plot(these_last_eig_below_thr(ii),these_angles(these_last_eig_below_thr(ii),ii),'marker','.','markersize',38,...
            'color',colors(i,:),'linestyle','none')
    end
end
set(gca,'Tickdir','out'),set(gca,'FontSize',14)
xlabel('dimension'),ylabel('principal angle (deg)')
xlim([0 dim_manifold+1]),ylim([0 90])



% -------------------------------------------------------------------------
% Histogram with the number of non-orthogonal pairs of eigenvectors per
% type of task

% matrix for stack plot
matrix_stack_hist_princ_angles = cell2mat( arrayfun( @(x) x.hist.counts_below_th, ...
                                summary_data.princ_angles, 'UniformOutput', false )' );
                            
% legend
legend_plot                 = cell(numel(meta_info.task_pairs.unique_pairs),1);
for i = 1:numel(meta_info.task_pairs.unique_pairs)
    legend_plot{i}          = [meta_info.task_pairs.unique_pairs{i}{1} ' vs. ' ...
                                meta_info.task_pairs.unique_pairs{i}{2}];
end

figure,
bar(summary_data.princ_angles(1).hist.edges', matrix_stack_hist_princ_angles','stack','EdgeColor',[1 1 1])
set(gca,'TickDir','out'), set(gca,'FontSize',16)
ylim([0 max(sum(matrix_stack_hist_princ_angles))+1]), xlim([0 dim_manifold+1])
ylabel('counts'),xlabel(['number non-orthogonal dimensions (P<' num2str(Pval_empir_angle_dist) ')'])
legend(legend_plot,'Location','NorthWest'), legend boxoff