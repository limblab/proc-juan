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
dims = 1:20;

% Number of manifold dimensions for scatter plot that compares weights
% across all task pairs. Only for visualization
dims_scatter = 1:20;

% Dataset type
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFINE WHAT IS A UNIT WITH HIGH WEIGHT: Take the X percentile of PCA
% weights across all datasets

all_w = [];
all_wrist_u = [];

for d = 1:length(datasets)
    
     % retrieve PCA weights
     for t = 1:length(datasets{d}.dim_red_FR)

         if ~isempty(dims)
             w{t} = datasets{d}.dim_red_FR{t}.w(:,dims);
         else
             w{t} = datasets{d}.dim_red_FR{t}.w;
         end
     end

     % get distribution of PCA weights
     all_w = [all_w; reshape(cell2mat(w),[],1)]; %#ok<*AGROW>
     
     % and whether it's a wrist or reach-to-grasp dataset
     if ismember(d,wrist_ds)
         all_wrist_u = [all_wrist_u; ones(numel(cell2mat(w)),1)]; 
     else
         all_wrist_u = [all_wrist_u; zeros(numel(cell2mat(w)),1)]; 
     end
     clear w;
end

% compute threshold for considering a weight "high"
w_th = prctile(abs(all_w),100-p_high_w*100);


% cast to bool, for quick indexing
all_wrist_u = boolean(all_wrist_u);



% -----------------------
% Histogram of all weights and high-weight threshold
step_hist = 0.01;
x_axis_all_w = 0:step_hist:1+step_hist;
hist_all_w = histcounts(all_w,x_axis_all_w);
figure,hold on
bar(x_axis_all_w(1:end-1),100*hist_all_w/numel(all_w),...
        'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
yl = ylim; plot([w_th, w_th], yl, 'r', 'linewidth',1.5);
text(double(w_th+.05),double(yl(2)-1),['P<' num2str(p_high_w)],'Fontsize',14)
set(gca,'TickDir','out','FontSize',14), box off
if isempty(dims)
    xlabel(['Weights all modes'])
else
    xlabel(['Weights modes 1:' num2str(dims)])
end
ylabel('Percentage (%)'),title('Histogram across all tasks')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get all pairs of weights for all task comparisons


all_w_pairs = [];
all_w_pairs_wrist_flg = [];

for d = 1:length(datasets)

    
    % get all tasks pairs (order matters for this!)
    n_tasks = length(datasets{d}.labels);
    all_combs = nchoosek(1:n_tasks,2);
    all_combs = [all_combs; nchoosek(n_tasks:-1:1,2)];

    for c = 1:length(all_combs)

        t1 = all_combs(c,1);
        t2 = all_combs(c,2);

        % retrieve weights
        if isempty(dims_scatter)
            w1 = datasets{d}.dim_red_FR{t1}.w;
            w2 = datasets{d}.dim_red_FR{t2}.w;
        else
            w1 = datasets{d}.dim_red_FR{t1}.w(:,1:dims_scatter);
            w2 = datasets{d}.dim_red_FR{t2}.w(:,1:dims_scatter);
        end
        
        % fill the matrices
        rw1 = reshape(w1,[],1);
        rw2 = reshape(w2,[],1);
        
        all_w_pairs = [all_w_pairs; rw1, rw2];
        if ismember(d,wrist_ds)
            all_w_pairs_wrist_flg = [all_w_pairs_wrist_flg; ones(length(rw1),1)];
        else
            all_w_pairs_wrist_flg = [all_w_pairs_wrist_flg; zeros(length(rw1),1)];
        end
    end
end

all_w_pairs_wrist_flg = boolean(all_w_pairs_wrist_flg);


% -----------------------
% Scatter plot of all pairs of weights for all modes in all task
% comparisons

col_wrist = [247 152 80]/255;
col_reach = [157 157 182]/255;

% This one doesn't take the absolute value
figure,hold on
plot(all_w_pairs(all_w_pairs_wrist_flg,1),all_w_pairs(all_w_pairs_wrist_flg,2),...
    '.','linestyle','none','color',col_wrist)
plot(all_w_pairs(~all_w_pairs_wrist_flg,1),all_w_pairs(~all_w_pairs_wrist_flg,2),...
    '.','linestyle','none','color',col_reach)
% high-weight threshold
plot([-1 -w_th],[w_th w_th],'r','linewidth',2)
plot([w_th 1],[w_th w_th],'r','linewidth',2)
text(-.9,.9,'High weight','color','r','FontSize',12)
plot([-1 -w_th],[-w_th -w_th],'r','linewidth',2)
plot([w_th 1],[-w_th -w_th],'r','linewidth',2)
plot([-w_th -w_th],[-1 -w_th],'r','linewidth',2)
plot([-w_th -w_th],[w_th 1],'r','linewidth',2)
plot([w_th w_th],[-1 -w_th],'r','linewidth',2)
plot([w_th w_th],[w_th 1],'r','linewidth',2)
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
xlim([-1 1]),ylim([-1 1])
xlabel('Weight task 1'),ylabel('Weight task 2')
if isempty(dims_scatter)
    title('PCA weights across tasks for all neural modes')
else
    title(['PCA weights across tasks for all the leading ' num2str(dims_scatter(end)) ' modes'])
end


% And this one does take the absolute value
figure,hold on
plot(abs(all_w_pairs(all_w_pairs_wrist_flg,1)),abs(all_w_pairs(all_w_pairs_wrist_flg,2)),...
    '.','linestyle','none','color',col_wrist)
plot(abs(all_w_pairs(~all_w_pairs_wrist_flg,1)),abs(all_w_pairs(~all_w_pairs_wrist_flg,2)),...
    '.','linestyle','none','color',col_reach)
% high-weight threshold
plot([w_th w_th],[w_th 1],'r','linewidth',2)
plot([w_th 1],[w_th w_th],'r','linewidth',2)
text(.8,double(w_th+.05),'High weight','color','r','FontSize',12)
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
xlim([0 1]),ylim([0 1])
xlabel('Weight task 1'),ylabel('Weight task 2')
if isempty(dims_scatter)
    title('PCA weights across tasks for all neural modes')
else
    title(['PCA weights across tasks for all the leading ' num2str(dims_scatter(end)) ' modes'])
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For each pair of tasks, CHECK WHAT UNITS HAVE BOTH WEIGHTS IN BOTH
% -- For this analysis we can't do all pair-wise comparisons, because the
% order does matter: a unit can have a high weight in task i but not in
% task j and viceversa


ctr = 1;

for d = 1:length(datasets)

    
    % get all tasks pairs (order matters for this!)
    n_tasks = length(datasets{d}.labels);
    all_combs = nchoosek(1:n_tasks,2);
    all_combs = [all_combs; nchoosek(n_tasks:-1:1,2)];

    for c = 1:length(all_combs)

        t1 = all_combs(c,1);
        t2 = all_combs(c,2);

        % retrieve weights
        if isempty(dims)
            w1 = datasets{d}.dim_red_FR{t1}.w;
            w2 = datasets{d}.dim_red_FR{t2}.w;
        else
            w1 = datasets{d}.dim_red_FR{t1}.w(:,1:dims);
            w2 = datasets{d}.dim_red_FR{t2}.w(:,1:dims);
        end

        % ---------------------
        % Find units with high weights in task 1
        [ro, co] = find(w1>w_th);
        [uniq_ro, idx_uniq_ro] = unique(ro);

        % ---------------------
        % See if these same units have a high weight contribution to any
        % mode in task 2

        [ro2, co2] = find(w2>w_th);
        [uniq_ro2, idx_uniq_ro2] = unique(ro2);

        % high-weight units in both tasks
        if length(uniq_ro)<length(uniq_ro2)
            u_high_w_both = uniq_ro(ismember(uniq_ro,uniq_ro2));
        else
            u_high_w_both = uniq_ro2(ismember(uniq_ro2,uniq_ro));
        end

    %     % figure to check that things look right
    %     figure,
    %     subplot(121),imagesc(w1>w_th),hold on
    %     plot([ones(length(u_high_w_both),1)*1, ones(length(u_high_w_both),1)*size(w1,2)]',[u_high_w_both, u_high_w_both]','r')
    %     subplot(122),imagesc(w2>w_th),hold on
    %     plot([ones(length(u_high_w_both),1)*1, ones(length(u_high_w_both),1)*size(w1,2)]',[u_high_w_both, u_high_w_both]','r')


        % for each unit with a high-weight in both tasks, get the modes that
        % have high weight in each task
        [r1, c1] = find(w1(u_high_w_both,:)>w_th);
        [r2, c2] = find(w2(u_high_w_both,:)>w_th);


        % save the modes in task 1 that have high weight for each unit and
        % their weight 

        modes_high_w1 = [];
        w_modes_high_w1 = [];
        modes_high_w2 = [];
        w_modes_high_w2 = [];
        
        for j = 1:length(u_high_w_both)
            modes_high_w1{j} = find(w1(u_high_w_both(j),:)>w_th);
            w_modes_high_w1{j} = w1(u_high_w_both(j),modes_high_w1{j});

            modes_high_w2{j} = find(w2(u_high_w_both(j),:)>w_th);
            w_modes_high_w2{j} = w2(u_high_w_both(j),modes_high_w2{j});
        end


        % Store results
        all_high_w(ctr).high_w_units = u_high_w_both; %#ok<*SAGROW>
        all_high_w(ctr).n_high_w_units = length(u_high_w_both); %#ok<*SAGROW>
        all_high_w(ctr).modes1 = w1(u_high_w_both,:);
        all_high_w(ctr).modes2 = w2(u_high_w_both,:);
        all_high_w(ctr).modes_high_w1 = modes_high_w1;
        all_high_w(ctr).modes_high_w2 = modes_high_w2;
        all_high_w(ctr).w_modes_high_w1 = w_modes_high_w1;
        all_high_w(ctr).w_modes_high_w2 = w_modes_high_w2;
        all_high_w(ctr).ds = d;
        all_high_w(ctr).comp = [t1 t2];
        all_high_w(ctr).n_units = size(w1,1);
        if ismember(d,wrist_ds)
            all_high_w(ctr).wrist_ds = true;
        else
            all_high_w(ctr).wrist_ds = false;
        end

        ctr = ctr + 1;
    end
end