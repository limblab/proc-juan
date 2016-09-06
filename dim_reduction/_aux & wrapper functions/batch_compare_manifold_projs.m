%
%
%

function proj_results = batch_compare_manifold_projs( varargin )


% -------------------------------------------------------------------------
% read input params

if nargin >= 1
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end
if nargin == 2
    params              = varargin{2};
else
    params              = batch_compare_manifold_projs_defaults();
end


% -------------------------------------------------------------------------
% load data

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end


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


% -------------------------------------------------------------------------
% get info about monkeys, tasks, and who did what

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% compare manifold projections using CCA

data                    = cell(1,length(meta_info.tasks_per_session));

for i = 1:meta_info.nbr_monkeys
    
    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
        
        % what dataset are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);

        
        % -----------------------------
        % compare the projections by doing CCA between all pairs of tasks
        % in this session 
        can_corrs       = canon_corr_all_manifolds( datasets{dtst}.stdata, ...
                            1:params.dim_manifold, datasets{dtst}.labels, ...
                            params.target, params.time_win(dtst,:) );
                        

        % -----------------------------
        % store all results
        data{dtst}.can_corrs = can_corrs;
    end
end




% -------------------------------------------------------------------------
% Do some summary analyses

% how many pairs of combinations ( monkey x session x pairs of tasks ) do we have?
nbr_manifold_pairs      = length(meta_info.task_pairs.monkey);
nbr_pairs_tasks      	= length(meta_info.task_pairs.unique_pairs);


% create struct to store the summary results
summary_data.canon_corrs = struct('pair',[meta_info.task_pairs.unique_pairs],...
                            'nbr_canon_corrs_above_chance',[],...
                            'cc',zeros(1,params.dim_manifold));

% create array with all canonical correlations
cc_matrix               = zeros(nbr_manifold_pairs,params.dim_manifold);

for d = 1:length(can_corrs)
    for p = 1:size(data{d}.can_corrs.cc,1)
        
    end
end
                        
                        
% ------------------------
% summarize, for each pair of tasks (e.g., tube vs ball, or iso vs mov) the
% canonical correlation for each pair of projections
for i = 1:nbr_manifold_pairs
   for ii = 1:nbr_pairs_tasks
        % retrieve pairs of tasks
        pair_cmp        = strcmpi(meta_info.task_pairs.unique_pairs{ii},...
                            sort(meta_info.task_pairs.task_pair{i}));

        if sum(pair_cmp) == 2
            % store the canoncorrs per task. Just append each trial one
            % after the other
            if numel(unique(summary_data.canon_corrs(ii).cc)) == 1 
                % summary_data.canon_corrs(ii).cc(1,:) = 
            else
            end
        end
   end
end



% -------------------------------------------------------------------------
% return variable
proj_results.data = can_corrs;
proj_results.summary_data = summary_data;


% -------------------------------------------------------------------------
% PLOTS

% ---------------------------------------------
% canonical correlations for all pairs of tasks
colors                      = parula(nbr_pairs_tasks);

figure, hold on
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