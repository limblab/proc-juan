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

function batch_dimensionality_analysis( varargin )


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

meta_info = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% summarize variance explained by different means


colors_p_task           = parula(meta_info.nbr_tasks);


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


% one figure with all meta_info.monkeys and meta_info.tasks
linestyle_p_monkey = {'-',':','-.','--'};

figure, hold on
for i = 1:meta_info.nbr_monkeys
    for ii = 1:numel(meta_info.sessions_per_monkey{i})
        for iii = 1:length(meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)})
            % see what task it is
            this_task   = meta_info.tasks_per_session{meta_info.sessions_per_monkey{i}(ii)}(iii);
            % and assign color
            this_color  = colors_p_task(find(strncmpi(meta_info.tasks, this_task, length(this_task) ),1),:);
            % plot
            plot(cumsum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)/...
                sum(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen),...
                'color',this_color,'linewidth',2,'linestyle',linestyle_p_monkey{i})
            set(gca,'TickDir','out'), ylim([0 1]),set(gca,'FontSize',16)
%            xlim([0 length(datasets{meta_info.sessions_per_monkey{i}(ii)}.dim_red_FR{iii}.eigen)]);
            xlim([0 50])
            xlabel('dimension'), ylabel('variance explained');
        end
    end
end
