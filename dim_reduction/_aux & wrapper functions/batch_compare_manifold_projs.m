%
%
%

function proj_results = batch_compare_manifold_projs( varargin )


% check that the parallel pool is running, otherwise start it
gcp;


% -------------------------------------------------------------------------
% read inputs

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
    params              = batch_compare_manifold_projs_defaults( varargin{2} );
else
    params              = batch_compare_manifold_projs_defaults();
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

        disp(['comparing neural projs. dataset #' num2str(dtst)]);

        
        % -----------------------------
        % compare the projections by doing CCA between all pairs of tasks
        % in this session 
        can_corrs       = canon_corr_all_manifolds( datasets{dtst}.stdata, ...
                            1:params.dim_manifold, datasets{dtst}.labels, ...
                            params.target, params.time_win(dtst,:) );
                        
        
        % get within task CC by taking 100 random subsets of half of the
        % trials and computing the CC across them
        cc_within       = canon_corr_within_all_manifolds( datasets{dtst}.stdata, ...
                            1:params.dim_manifold, true, 'all_conc', ...
                            params.time_win(dtst,:), 100 );
                        
                        
        % -----------------------------
        % If want to obtain confidence limit with bootstrapping
        if params.do_bootstrap
            signif_boots = zeros(length(can_corrs.lin_transform),params.dim_manifold);
            for p = 1:length(can_corrs.lin_transform)
                % Bootstrapping shuffling in time
                signif_boots(p,:) = bootstrap_canon_corr( can_corrs.lin_transform(p).U, ...
                            can_corrs.lin_transform(p).V, params.nbr_shuffles_bootstrap,...
                            params.prctile_bootstrap );
            end
%             % Bootstrapping shuffling weights of neural units onto
%             % modes
%             signif_boots2 = bootstrap_weights_canon_corr( datasets{dtst}.stdata, ...
%                                 1:params.dim_manifold, 'all_conc', params.time_win(dtst,:),...
%                                 params.nbr_shuffles_bootstrap, params.prctile_bootstrap);

        end
                        
                        
        % -----------------------------
        % store all results
        
        % the CCA results
        data{dtst}.can_corrs = can_corrs;
        % the within task CC
        data{dtst}.within_can_corrs = cc_within;
        % bootstrapping results, if done
        if params.do_bootstrap
            data{dtst}.can_corrs.signif_boots = signif_boots;
%             % TO DELETE
%             data{dtst}.can_corrs.signif_boots_weights = signif_boots2;
            data{dtst}.can_corrs.params_boots.prctile_signif = params.prctile_bootstrap;
            data{dtst}.can_corrs.params_boots.nbr_shuffles = params.nbr_shuffles_bootstrap;
        end
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
                            'cc',zeros(1,params.dim_manifold),...
                            'norm_cc',zeros(1,params.dim_manifold),...
                            'ceil_cc',zeros(1,params.dim_manifold));


                        
% ------------------------
% Get P value obtained with Matlab's built-in stat test; in case haven't
% decided not to do bootstrapping
pooled_stats            = cell2mat(cellfun(@(y) y.can_corrs.stats, data,...
                                'UniformOutput',false));
pooled_P_vals           = cell2mat(arrayfun(@(x) x.p, pooled_stats,'UniformOutput',false)');


% ------------------------
% Fill 2D matrices with results, so it's easier to compile them

% matrix with all canonical correlations
nbr_projs_above_chance  = zeros(1,nbr_manifold_pairs);
cc_matrix               = zeros(nbr_manifold_pairs,params.dim_manifold);
chance_matrix           = zeros(nbr_manifold_pairs,params.dim_manifold);
norm_cc_matrix          = zeros(nbr_manifold_pairs,params.dim_manifold);
ceil_cc_matrix          = zeros(nbr_manifold_pairs,params.dim_manifold);
ctr                     = 1;
    
for d = 1:length(data)
    for p = 1:size(data{d}.can_corrs.cc,1)
        
        % fill matrix with all CCs
        cc_matrix(ctr,:) = data{d}.can_corrs.cc(p,:);
        
        % fill chance levels
        chance_matrix(ctr,:) = data{d}.can_corrs.signif_boots(p,:);
        
        % ceil CC (max of the 99th percentiel of the within CC)
        % -- it used to be the mean but it was changed following a
        % reviewer's comment
        ceil_cc_matrix(ctr,:) = data{d}.within_can_corrs.pair{p}.pctile99_max;
        
        % norm CC matrix
        norm_cc_matrix(ctr,:) = cc_matrix(ctr,:)./chance_matrix(ctr,:);
        
        % get number of projections that are greater than change
        if params.do_bootstrap
            nbr_projs_above_chance(ctr) = sum( cc_matrix(ctr,:) > chance_matrix(ctr,:) );
        else
            nbr_projs_above_chance(ctr) = sum( pooled_P_vals(ctr,:) < params.P_thr );
        end
        
        ctr             = ctr + 1;
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
                summary_data.canon_corrs(ii).cc(1,:) = cc_matrix(i,:);
            else
                summary_data.canon_corrs(ii).cc = [summary_data.canon_corrs(ii).cc;
                    cc_matrix(i,:)];
            end
            % store normalized CCs
            if numel(unique(summary_data.canon_corrs(ii).norm_cc)) == 1 
                summary_data.canon_corrs(ii).norm_cc(1,:) = norm_cc_matrix(i,:);
            else
                summary_data.canon_corrs(ii).norm_cc = [summary_data.canon_corrs(ii).norm_cc;
                    norm_cc_matrix(i,:)];
            end
            % store ceiling CCs (mean percentile99 across both tasks in the pair)
            if numel(unique(summary_data.canon_corrs(ii).ceil_cc)) == 1 
                summary_data.canon_corrs(ii).ceil_cc(1,:) = ceil_cc_matrix(i,:);
            else
                summary_data.canon_corrs(ii).ceil_cc = [summary_data.canon_corrs(ii).ceil_cc;
                    ceil_cc_matrix(i,:)];
            end
            % get nbr projections whose correlations are above chance
            summary_data.canon_corrs(ii).nbr_canon_corrs_above_chance = [...
                summary_data.canon_corrs(ii).nbr_canon_corrs_above_chance, ...
            nbr_projs_above_chance(i)];
        
            break;
        end
   end
end



% -------------------------------------------------------------------------
% return variable
proj_results.data           = data;
proj_results.summary_data   = summary_data;
proj_results.params         = params;


% -------------------------------------------------------------------------
% PLOTS

% ---------------------------------------------
% canonical correlations for all pairs of tasks
colors                      = parula(nbr_pairs_tasks);

figure, hold on
for i = 1:nbr_pairs_tasks
    these_ccs               = summary_data.canon_corrs(i).cc';
    % what projections are above the chance level according to matlab's
    % test?
    these_proj_below_chane  = summary_data.canon_corrs(i).nbr_canon_corrs_above_chance;
    
    % plot the traces
	plot(these_ccs,'color',colors(i,:),'linewidth',2)
    % and add markers showing the last dim below the randomness threshold
    for ii = 1:size(these_ccs,2)
        plot(these_proj_below_chane(ii),these_ccs(these_proj_below_chane(ii),ii),'marker','.','markersize',38,...
            'color',colors(i,:),'linestyle','none')
    end
end
set(gca,'Tickdir','out'),set(gca,'FontSize',14)
xlabel('projection'),ylabel('canonical correlation')
xlim([0 params.dim_manifold+1]),ylim([0 1])


% ---------------------------------------------
% Normalized Canonical correlation (divided by the significance level)

figure, hold on
for i = 1:nbr_pairs_tasks
    these_norm_ccs          = summary_data.canon_corrs(i).norm_cc';
    
    % plot the traces
	plot(these_norm_ccs,'color',colors(i,:),'linewidth',2)
end
plot([1 params.dim_manifold],[1 1],'color',[.6 .6 .6],'linewidth',2,'linestyle','-.')
set(gca,'Tickdir','out'),set(gca,'FontSize',14)
xlabel('projection'),ylabel('normalized canonical correlation')
xlim([0 params.dim_manifold+1]),ylim([0 ceil(max(max(norm_cc_matrix)))])



% ---------------------------------------------
% 1 plot with canonical correlations for all the task the monkey performed
% in that session (i.e.) one plot per task

if params.plot_p_session
    for i = 1:length(meta_info.tasks_per_session)

        these_ccs               = data{i}.can_corrs.cc';
        these_ceilings          = cellfun( @(x) x.pctile99_max, data{i}.within_can_corrs.pair, ...
                                    'UniformOutput', false );
        these_ceilings          = cell2mat( these_ceilings' )';
        
        % get bootstrapping traces if they were computed, otherwise use
        % matlab's canon_corr significance test
        if isfield(data{1}.can_corrs,'signif_boots') 
            % read bootstrap correlations
            these_bootstrap     = data{i}.can_corrs.signif_boots';
        else
            % find what is the last projection whose correlation is higher than
            % chance
            these_P_chance      = cell2mat(arrayfun( @(x) x.p, data{i}.can_corrs.stats,...
                                    'UniformOutput',false )');
            these_above_P_chance =  these_P_chance > params.P_thr;
            these_last_proj_below_P_chance = repmat(params.dim_manifold,size(these_ccs,2),1) ...
                    - sum(these_above_P_chance,2) - ones(size(these_ccs,2),1);
        end

        % get what task pair these are, in the big scheme of things, to choose
        % the color for plotting and make the legend
        these_task_pairs        = meta_info.task_pairs.task_pair_nbr( meta_info.task_pairs.session == i );
        % make legend
        for p = 1:length(these_task_pairs)
            this_legend{p}         = [meta_info.task_pairs.task_pair{these_task_pairs(p)}{1} ' vs ' ...
                                        meta_info.task_pairs.task_pair{these_task_pairs(p)}{2}];
        end
        if exist('these_bootstrap','var')
            this_legend{length(this_legend)+1} = 'bootstrapped';
        end

        % plot!
        figure, hold on
        for p = 1:size(these_ccs,2)
            plot(these_ccs(:,p),'linewidth',2,'color',colors(these_task_pairs(p),:))
            if ~exist('these_bootstrap','var')
                plot(these_last_proj_below_P_chance(p),these_ccs(these_last_proj_below_P_chance(p),p),...
                    'marker','.','markersize',24,'color',colors(these_task_pairs(p),:))
            end
        end
        if exist('these_bootstrap','var')
            plot(these_bootstrap(:,p),':','linewidth',2,'color',[.6 .6 .6])
        end
        legend(this_legend,'Location','NorthEast'), legend boxoff
        set(gca,'Tickdir','out'),set(gca,'FontSize',14)
        xlim([0 params.dim_manifold+1]),ylim([0 1])
        xlabel('projection'),ylabel('canonical correlation')

        % plot ceilings
        for p = 1:size(these_ceilings,2)
            plot(these_ceilings(:,p),':','linewidth',2,'color',colors(these_task_pairs(p),:))
        end
        
        clear this_legend;
    end
end


% ---------------------------------------------
% Histogram with nbr of canonical correlations above chance

hist_matrix = cell2mat( arrayfun( @(x) histcounts(x.nbr_canon_corrs_above_chance,1:21), ...
                proj_results.summary_data.canon_corrs, 'UniformOutput', false )' );

figure, bar(1:20,hist_matrix','stacked')
set(gca,'TickDir','out','FontSize',14), xlim([0 params.dim_manifold+2])
xlabel('number signif. canon. corrs. (P<0.01)'),ylabel('counts')
% legend(meta_info.task_pairs.unique_pairs{1})


% 
% % ---------------------------------------------
% % P vals of CC
% 
% figure,plot(pooled_P_vals','color',[.6 .6 .6])
% set(gca,'TickDir','out'), set(gca,'FontSize',14);
% xlabel('projection'), ylabel('P val')
% ylim([0 .01])