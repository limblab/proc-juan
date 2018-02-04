% 
% Test if units with high PCA weights on one task also have high weights on
% another
%
% It generates FIG. 3 "ARE UNITS WITH HIGH PCA WEIGHTS THE SAME ACROSS
% TASKS?" of the Nature Comms review letter, which was later incorporated as
% a panel in a Suppl Fig
%


% What percentile of the PCA weights we are going to consider as high
p_high_w = 0.01; % 1; % (all units) %0.01; % 0.05; 

% Number of manifold dimensions ---all (full space) if empty
n_dims = [];

% Dataset type
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];


% -------------------------------------------------------------------------
% Define what is a unit with high weight: Take the X percentile of PCA
% weights across all datasets

all_w = [];

for d = 1:length(datasets)
    
     % retrieve PCA weights
     for t = 1:length(datasets{d}.dim_red_FR)

         if ~isempty(n_dims)
             w{t} = datasets{d}.dim_red_FR{t}.w(:,1:n_dims);
         else
             w{t} = datasets{d}.dim_red_FR{t}.w;
         end
     end

     % get distribution of PCA weights
     all_w = [all_w; reshape(cell2mat(w),[],1)];
     
     clear w;
end

% compute threshold for considering a weight "high"
w_th = prctile(abs(all_w),100-p_high_w*100);


% -------------------------------------------------------------------------


% to store all pairs of units in which the first has high weight
all_high_w_pairs = [];
% bool that tells if the pair above is from a wrist (1) or a reach-to-grasp
% (0) task 
wrist_u = [];


% What dataset?
for d = 1:length(datasets)

    % retrieve PCA weights for all the tasks in this dataset
    for t = 1:length(datasets{d}.dim_red_FR)

        if ~isempty(n_dims)
            w{t} = datasets{d}.dim_red_FR{t}.w(:,1:n_dims);
        else
            w{t} = datasets{d}.dim_red_FR{t}.w;
        end
    end

    
%     % ---------------------------------------------------------------------
%     % THIS IS WEIRD - I THINK IT'S WRONG !!!! --CHECK
%     
%     % get distribution of PCA weights
%     all_w = reshape(cell2mat(w),[],1);
% 
%     % compute threshold for considering a weight "high"
%     w_th = prctile(abs(all_w),100-p_high_w*100);


    % -----------------------
    % Find the units above that threshold, and whether they are the same across
    % different tasks

    % get all combinations of tasks
    comb_t = nchoosek(1:length(w),2);

    for c = 1:size(comb_t,1)

        % find indexes of units with high weights --the row (ro) is the unit
        % number
        [ro, co] = find(w{comb_t(c,1)}>=w_th);

        for u = 1:length(ro)
            w1 = abs(w{comb_t(c,1)}(ro(u),co(u)));
            w2 = abs(w{comb_t(c,2)}(ro(u),co(u)));
            
            % store whether it's a wrist or reach-to-grasp dataset
            if ismember(d,wrist_ds), wrist_u = [wrist_u, 1]; else 
                wrist_u = [wrist_u, 0]; end
            
            % store pair of weights
            all_high_w_pairs = [all_high_w_pairs; w1, w2];
        end
    end

    clearvars -except datasets n_dims p_high_w all_high_w_pairs wrist_u *_ds w_th
end

% make boolean, for quick indexing
wrist_u = boolean(wrist_u);


% -----------------------
% Scatter plot weights of high-weight units in task 1 vs weights in task 2

figure,hold on
plot(all_high_w_pairs(wrist_u,1),all_high_w_pairs(wrist_u,2),'.k','markersize',12)
plot(all_high_w_pairs(~wrist_u,1),all_high_w_pairs(~wrist_u,2),'.','color',[.6 .6 .6],'markersize',12)
xlim([0 1]),ylim([0 1])
xlabel('Task 1'),ylabel('Task 2')
set(gca,'TickDir','out','FontSize',14)
title(['Cross-task abs(PCA weights) >' num2str(100-p_high_w*100) '% percentile'])



% -----------------------
% Histograms with corr weights per task

x_hist = 0:0.05:1;
y_wrist = histcounts(all_high_w_pairs(wrist_u,2),x_hist)/numel(find(wrist_u));
y_reach = histcounts(all_high_w_pairs(~wrist_u,2),x_hist)/numel(find(~wrist_u));

% compute how many units with high weights also have high weights in
% another task
perc_wrist = numel(find(all_high_w_pairs(wrist_u,2)>w_th))/numel(all_high_w_pairs(wrist_u,2)>w_th);
perc_reach = numel(find(all_high_w_pairs(~wrist_u,2)>w_th))/numel(all_high_w_pairs(~wrist_u,2)>w_th);

figure,
subplot(121), hold on, 
bar(x_hist(1:end-1),y_wrist*100,'FaceColor',[0 0 0]), title(['Wrist datasets']),
plot([w_th w_th],[0 max(y_wrist)*100],'r','linewidth',2)
set(gca,'TickDir','out','FontSize',14), xlabel('Correlation'), ylabel('Percentage (%)')
text(0.4,floor(max(y_wrist)*10)*10,['n=' num2str(numel(find(wrist_u)))],'Fontsize',14);
text(0.4,floor(max(y_wrist)*10)*10-5,['% high-w=' num2str(perc_wrist*100,3)],'Fontsize',14);
subplot(122), hold on
bar(x_hist(1:end-1),y_reach*100,'FaceColor',[.6 .6 .6]), title('Reach-to-grasp datasets'), 
plot([w_th w_th],[0 max(y_reach)*100],'r','linewidth',2)
set(gca,'TickDir','out','FontSize',14), xlabel('Correlation')
text(0.4,floor(max(y_reach)*10)*10,['n=' num2str(numel(find(~wrist_u)))],'Fontsize',14);
text(0.4,floor(max(y_reach)*10)*10-5,['% high-w=' num2str(perc_reach*100,3)],'Fontsize',14);

