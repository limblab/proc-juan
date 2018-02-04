%
% Check dimensionality per task type. 
%
% BIG ISSUE: Something is not working right when doing PCA of the Ball task
% -- is it because the trial-averaged matrix is rank defficient ????
%


ds_to_use = [1:3 7:9];


% compute the dimensionality

ctr = 1;

for ds = 1:length(ds_to_use)
    
    for t = 1:length(datasets{ds_to_use(ds)}.labels)
        
        nfa = noise_floor_pca(datasets{ds_to_use(ds)}.stdata{t});
        
        % store dimensionality
        dimensionality(ctr) = nfa.dims;
        % and task and monkey
        task(ctr).label = datasets{ds_to_use(ds)}.labels{t};
        monkey(ctr).name = datasets{ds_to_use(ds)}.monkey;
        
        ctr = ctr + 1;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the distribution of task dimensionality


tasks = unique({task.label});
x_axis = 1:11;
h_dim = [];

for t = 1:length(tasks)
    
    % find instances of this task
    t_dim = dimensionality(strcmp({task.label},tasks(t)));
    
    % histogram
    h_dim(t,:) = histcounts(t_dim,x_axis);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT

% Dimensionality per task
figure
bar(x_axis(1:end-1),h_dim','Stacked')
set(gca,'TickDir','out'),set(gca,'FontSize',14); box off
ylabel('Counts'),xlabel('Estimated manifold dimensionality')
legend(tasks,'Location','NorthEast'), legend boxoff
yl = ylim;
text(1, yl(2)-.5, ['n = ' num2str(length(dimensionality))],'FontSize',14)


% Dimensionality for all tasks combined
figure
bar(x_axis(1:end-1),sum(h_dim,1),'FaceColor',[.5 .5 .5])
set(gca,'TickDir','out'),set(gca,'FontSize',14); box off
ylabel('Counts'),xlabel('Estimated manifold dimensionality')
yl = ylim;
text(1, yl(2)-.5, ['n = ' num2str(length(dimensionality))],'FontSize',14)
