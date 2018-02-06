%
% Understanding the "dimensionality" of the different tasks by looking at the
% eigenvalue distribution using some simple ideas
%


% Params
nbr_iter = 100;


% get total number of tasks
nbr_tasks = sum( cellfun(@(x) length(x.labels), datasets) );

% Find minimum number of units across all datasets
nbr_chs = cellfun(@(x) length(x.neural_chs), datasets);
min_nbr_chs = min(nbr_chs);


%% ------------------------------------------------------------------------
% Take random subsets of 'min_nbr_chs', do PCA and store the mean and SD
% eigenvalue distribution 

m_eigenv = zeros(nbr_tasks,min_nbr_chs);
sd_eigenv = zeros(nbr_tasks,min_nbr_chs);
m_cumsum_eigenv = zeros(nbr_tasks,min_nbr_chs);
sd_cumsum_eigenv = zeros(nbr_tasks,min_nbr_chs);
ctr = 1;

monkey = [];
task = [];

% do for each dataset
for d = 1:length(datasets)
    
    
    % do for each task
    for t = 1:length(datasets{d}.labels)
        
        % get FRs of all neurons
        fr = datasets{d}.stdata{t}.target{end}.neural_data.conc_smoothed_fr;
        
        t_eigenv = zeros(nbr_iter,min_nbr_chs);
        
        % ----------------------
        % do 'nbr_iter' times 
        for i = 1:nbr_iter
            
            % take a subset of channels
            t_chs = datasample(1:length(datasets{d}.neural_chs),min_nbr_chs,'Replace',false);
    
            [~,~,t_eigenv(i,:)] = pca(fr(:,t_chs));
        end
        
        % take mean and SD of the eigenvalue distributions
        m_eigenv(ctr,:) = mean(t_eigenv,1);
        sd_eigenv(ctr,:) = std(t_eigenv,0,1);
        m_cumsum_eigenv(ctr,:) = mean(cumsum(t_eigenv,2)./sum(t_eigenv,2));
        sd_cumsum_eigenv(ctr,:) = std(cumsum(t_eigenv,2)./sum(t_eigenv,2),0,1);
        
        monkey{ctr} = datasets{d}.monkey;
        task{ctr} = datasets{d}.labels{t};
        
        ctr = ctr + 1;
    end
end
        
        
%% ------------------------------------------------------------------------

% SOME BASIC PLOTS


col.ball = [128 0 128]/255; % purple
col.mg_pt = [157 157 182]/255; % violet
col.iso = [255 0 0]/255; % red
col.wm = [255,165,0]/255; % orange
col.spr = [255,215,0]/255; % gold
col.iso8 = [128,0,0]/255; % maroon (dark red)

col.C_wrist  = [255 165 0]/255;
col.J_wrist = [1 0 0];
col.C_reach = [.6 .6 .6];
col.T_reach = [0 0 0];


% Need to rename all mg-pt to mg_pt so the code below works...
task(strcmp(task,'mg-pt')) = {'mg_pt'};

% -------------------------------------------------------------------------
% CUMSUM OF THE VAF PER MONKEY AND DATASET TYPE
% Need to improve this code

figure, hold on
plot(m_cumsum_eigenv(1,:)','color',col.C_wrist,'linewidth',1.5)
plot(m_cumsum_eigenv(16,:)','color',col.J_wrist,'linewidth',1.5)
plot(m_cumsum_eigenv(10,:)','color',col.C_reach,'linewidth',1.5)
plot(m_cumsum_eigenv(27,:)','color',col.T_reach,'linewidth',1.5)

plot(m_cumsum_eigenv(2:9,:)','color',col.C_wrist,'linewidth',1.5)
plot(m_cumsum_eigenv(17:26,:)','color',col.J_wrist,'linewidth',1.5)
plot(m_cumsum_eigenv(11:15,:)','color',col.C_reach,'linewidth',1.5)
plot(m_cumsum_eigenv(28:30,:)','color',col.T_reach,'linewidth',1.5)

set(gca,'TickDir','out','FontSize',14), ylim([0 1])
xlabel('Neural mode'), ylabel ('Percentage variance explained (%)')
legend('C - wrist','J - wrist', 'C - reach','T - reach','Location','Southeast') 
legend boxoff, box off


% -------------------------------------------------------------------------
% % Plot SD eigenvalue distribution

% figure, hold on
% plot(100*sd_eigenv(1,:)./sum(m_eigenv(1,:))','color',[255 165 0]/255,'linewidth',1.5)
% plot(100*sd_eigenv(16,:)./sum(m_eigenv(16,:))','color','r','linewidth',1.5)
% plot(100*sd_eigenv(10,:)./sum(m_eigenv(10,:))','color',[.6 .6 .6],'linewidth',1.5)
% plot(100*sd_eigenv(27,:)./sum(m_eigenv(27,:))','color','k','linewidth',1.5)
% 
% plot(100*(sd_eigenv(2:9,:)./sum(m_eigenv(2:9,:),2))','color',[255 165 0]/255,'linewidth',1.5)
% plot(100*(sd_eigenv(17:26,:)./sum(m_eigenv(17:26,:),2))','color','r','linewidth',1.5)
% plot(100*(sd_eigenv(11:15,:)./sum(m_eigenv(11:15,:),2))','color',[.6 .6 .6],'linewidth',1.5)
% plot(100*(sd_eigenv(28:30,:)./sum(m_eigenv(28:30,:),2))','color','k','linewidth',1.5)
% 
% set(gca,'TickDir','out','FontSize',14)
% xlabel('Neural mode'), ylabel ('SD across neural unit subsets - VAF (%)')
% legend('C - wrist','J - wrist', 'C - reach','T - reach','Location','Southeast') 
% legend boxoff, box off


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOTS PER MONKEY

u_monkeys = unique(monkey);


for m = 1:length(unique(monkey))

    % get what tasks these monkey performed
    i_task_t_monkey = find(strcmp(monkey,u_monkeys{m}));
    
    % get task names
    tasks_t_monkey = unique({task{i_task_t_monkey}});
    
    % ---------------------------------------------------------------------
    % plot cumsum eigenvals for each monkey
    % -- the two loops are to plot the legends correctly. There's a better
    % way, but I have no time for this now!
    figure, hold on
    for t = 1:length(tasks_t_monkey)
    
        % get indexes for this task
        i_wrt_task_t_monkey = find(strcmp({task{i_task_t_monkey}},tasks_t_monkey{t})); %#ok<CCAT1>
        i_abs = i_task_t_monkey(i_wrt_task_t_monkey); % need to do this to reference it to the matrix where everything is
        
        % get color from list in 'col'
        t_color = col.(tasks_t_monkey{t});

        % plot
        plot(100*m_cumsum_eigenv(i_abs(1),:)','linewidth',1.5,'color',t_color);
    end
    % this second loop is part of the dirty legend trick 
    for t = 1:length(tasks_t_monkey)
        
        % get indexes for this task
        i_wrt_task_t_monkey = find(strcmp({task{i_task_t_monkey}},tasks_t_monkey{t})); %#ok<CCAT1>
        i_abs = i_task_t_monkey(i_wrt_task_t_monkey); % need to do this to reference it to the matrix where everything is
        
        % get color from list in 'col'
        t_color = col.(tasks_t_monkey{t});
        
        % plot
        plot(100*m_cumsum_eigenv(i_abs(2:end),:)','linewidth',1.5,'color',t_color);
    end
    set(gca,'TickDir','out','FontSize',14)
    ylim([0 100]), title(['Dataset ' u_monkeys{m}]); box off
    xlabel('Neural mode'), ylabel ('Percentage variance explained (%)')
    legend(tasks_t_monkey,'Location','SouthEast','Interpreter','none'), legend boxoff
    xlim([0 12])
end
