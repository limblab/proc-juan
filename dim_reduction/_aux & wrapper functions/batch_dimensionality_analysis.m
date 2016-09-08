%
% Load the datasets in "all_manifold_datasets.mat", and plot the scree
% plots for each monkey and session. The code will also run Machen's method
% to estimate the variance due to neural noise and find an upper bound to
% the dimensionality of the task. --> ToDo
%
%   function batch_dimensionality_analysis( varargin )
%
%
% function batch_dimensionality_analysis( )
% function batch_dimensionality_analysis( all_datasets )
% function batch_dimensionality_analysis( path )
%
%

function dim_results = batch_dimensionality_analysis( varargin )


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
% get some meta data prior to calculating stuff & plotting

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% summarize variance explained in different ways

% mean explained variance per task for each monkey

dim_results         	= cell(1,meta_info.nbr_monkeys);

% for each monkey
for i = 1:meta_info.nbr_monkeys
    
    % S1) save monkey name
    dim_results{i}.monkey = meta_info.monkeys{i};
    dim_results{i}.tasks = cell(1,length(meta_info.tasks_per_monkey{i}));
    
    % and for each task it performed
    for ii = 1:length(meta_info.tasks_per_monkey{i})
        % get he task
        this_task       = meta_info.tasks_per_monkey{i}{ii};
        
        % S2) save this task
        dim_results{i}.tasks{ii} = this_task;
        
        % preallocate matrix to store explained variance --make it
        % arbitrarily to the number of sessions the monkey performed, and
        % 100 dimensions
        aux_eigenval    = zeros(100,meta_info.nbr_sessions_per_monkey(i));
        % look for its index in each session
        for iii = 1:meta_info.nbr_sessions_per_monkey(i)
            % get ptr to session
            this_session = meta_info.sessions_per_monkey{i}(iii);
            % get the indx of this task in this session
            task_indx   = find(strcmpi(meta_info.tasks_per_session{this_session},this_task));
            % now retrieve variance explained
            if ~isempty(task_indx)
                this_eigenval = datasets{this_session}.dim_red_FR{task_indx}.eigen;
                aux_eigenval(1:length(this_eigenval),iii) = this_eigenval';
            end
        end
        
        % delete empty row
        aux_eigenval(:,sum(aux_eigenval,1)==0) = [];
        
        % S4) save eigenvalues
        dim_results{i}.eigenval{ii} = aux_eigenval;
        % S5) save explained variance
        dim_results{i}.expl_var{ii} = cumsum(aux_eigenval)./repmat(sum(aux_eigenval),100,1);
    end
end




% -------------------------------------------------------------------------
% summarize the results

% find the number of dimensions that explain a % of the variance for each
% task and monkey

perc_var            = [.6 .7 .75 .8 .9];

for m = 1:meta_info.nbr_monkeys
    for t = 1:length(meta_info.tasks_per_monkey{m})        
        for r = 1:size(dim_results{m}.expl_var{t},2)
            for p = 1:length(perc_var)
                dim_results{m}.dims_var{t}(r,p) = find(dim_results{m}.expl_var{t}(:,r)>perc_var(p),1);
            end
        end
    end
    
    dim_results{m}.perc_var = perc_var;
end




% -------------------------------------------------------------------------
% Plots

% ---------------------------
% Number of dimensions to explain a % of the neural variance for each
% realization of each task

% get data for this percentage
perc_plot           = .75;
indx_perc_plot      = find(dim_results{1}.perc_var==perc_plot);
% results will be stored in a 2D matrix with nbr_tasks rows per monkey, and
% the monkeys ordered sequentially. The number of columns is 5, which is
% greater than the max nbr of repetitions of each task
aux_dim_perc        = zeros(meta_info.nbr_monkeys*meta_info.nbr_tasks,5);
aux_legend_ctr      = 1;

for m = 1:meta_info.nbr_monkeys
    for t = 1:length(meta_info.tasks)
        if sum(strcmpi(meta_info.tasks{t},meta_info.tasks_per_monkey{m})) == 1
            indx_task   = find(strcmpi(meta_info.tasks{t},dim_results{m}.tasks));
            if ~isempty(indx_task)
                aux_dims = dim_results{m}.dims_var{indx_task}(:,indx_perc_plot)';
                aux_dim_perc((m-1)*meta_info.nbr_tasks+t,1:length(aux_dims)) = aux_dims;
    %            aux_dim_perc
            end
            
            % create legend
            aux_legend{aux_legend_ctr} = meta_info.monkeys{m};
            aux_legend{aux_legend_ctr} = [aux_legend{aux_legend_ctr}, ' - ', meta_info.tasks{t}];
            aux_legend_ctr = aux_legend_ctr+1;
        end
    end
end

% remove rows of zeroes
aux_dim_perc(sum(aux_dim_perc,2)==0,:) = [];
% and cols of zeroes
aux_dim_perc(:,sum(aux_dim_perc,1)==0) = [];


% convert aux_dim_perc into histograms
x_dims_hist             = 1:20;
counts_dims_hist        = zeros(size(aux_dim_perc,1),length(x_dims_hist)-1);
for i = 1:size(aux_dim_perc,1)
    counts_dims_hist(i,:) = histcounts(aux_dim_perc(i,:),x_dims_hist);
end


figure,
bar(counts_dims_hist','stacked')
set(gca,'TickDir','out'), set(gca,'FontSize',16)
ylim([0 max(sum(counts_dims_hist,1)+1)])
xlim([0 max(x_dims_hist)])
xlabel(['dimension > ' num2str(perc_plot*100) ' % variance explained'])
ylabel('counts'), legend(aux_legend(1:end)), legend boxoff

% ---------------------------
% SCREE PLOTS --old code, very convoluted; can be easily improved

colors_p_task           = distinguishable_colors(meta_info.nbr_tasks,{'w'});


% plot each session that each monkey performed
for i = 1:meta_info.nbr_monkeys
    figure
    for ii = 1:numel(meta_info.sessions_per_monkey{i})
        subplot(1,numel(meta_info.sessions_per_monkey{i}),ii), hold on
        for iii = 1:length(meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)})
            % see what task it is
            this_task   = meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)}(iii);
            % and assign color
            this_color  = colors_p_task(find(strncmpi(meta_info.tasks, this_task, length(this_task) ),1),:);
            % plot
            plot(cumsum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)/...
                sum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen),...
                'color',this_color,'linewidth',2)
            set(gca,'TickDir','out'), ylim([0 1]),set(gca,'FontSize',12)
%             xlim([0 length(datasets{sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)]);
            xlim([0 20])
            xlabel('dimension')
        end
        legend(meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)},...
            'Location','SouthEast','FontSize',12)
        if ii == 1, 
            ylabel('variance explained'); 
            title([meta_info.monkeys{i}, ' ', datasets{meta_info.sessions_per_monkey{i}(ii)}.date(1:8), ...
                ' - n = ' num2str(length(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen))])  
        else
            title([datasets{meta_info.sessions_per_monkey{i}(ii)}.date(1:8), ' - n = ', ... 
                num2str(length(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen))]); 
        end
    end
end


% one figure per monkey, with all the tasks
for i = 1:meta_info.nbr_monkeys
    figure, hold on
    for ii = 1:numel(meta_info.sessions_per_monkey{i})
        for iii = 1:length(meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)})
            % see what task it is
            this_task   = meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)}(iii);
            % and assign color
            this_color  = colors_p_task(find(strncmpi(meta_info.tasks, this_task, length(this_task) ),1),:);
            % plot
            plot(cumsum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)/...
                sum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen),...
                'color',this_color,'linewidth',2)
            set(gca,'TickDir','out'), ylim([0 1]),set(gca,'FontSize',16)
            xlim([0 length(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)]);
%            xlim([0 20])
            xlabel('dimension'), title(meta_info.monkeys{i})
            legend(meta_info.tasks_per_monkey{i},'Location','SouthEast','FontSize',16)
            ylabel('variance explained');
        end
    end
end


% % one figure with all meta_info.monkeys and meta_info.tasks
% linestyle_p_monkey = {'-',':','-.','--'};
% 
% figure, hold on
% for i = 1:meta_info.nbr_monkeys
%     for ii = 1:numel(meta_info.sessions_per_monkey{i})
%         for iii = 1:length(meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)})
%             % see what task it is
%             this_task   = meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)}(iii);
%             % and assign color
%             this_color  = colors_p_task(find(strncmpi(meta_info.tasks, this_task, length(this_task) ),1),:);
%             % plot
%             plot(cumsum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)/...
%                 sum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen),...
%                 'color',this_color,'linewidth',2,'linestyle',linestyle_p_monkey{i})
%             set(gca,'TickDir','out'), ylim([0 1]),set(gca,'FontSize',16)
% %            xlim([0 length(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)]);
%             xlim([0 50])
%             xlabel('dimension'), ylabel('variance explained');
%         end
%     end
% end
