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

% Manifold dimensions (d_begin:d_end) ---all (full space) if empty. If
% negative (-X): last X dimensions
dims = [];

% Dimensions for finding the high-weight threshold
dims_th = [];

% Number of manifold dimensions for scatter plot that compares weights
% across all task pairs. Only for visualization
dims_scatter = [];

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

         if ~isempty(dims_th)
             if length(dims_th) > 1
                 w{t} = datasets{d}.dim_red_FR{t}.w(:,dims_th);
             else
                 w{t} = datasets{d}.dim_red_FR{t}.w(:,end+dims_th+1:end);
             end
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
hist_all_w = histcounts(abs(all_w),x_axis_all_w);

figure,hold on
bar(x_axis_all_w(1:end-1),100*hist_all_w/numel(all_w),...
        'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
set(gca,'Yscale','log')
yl = ylim; plot([w_th, w_th], yl, 'r', 'linewidth',1.5);
text(double(w_th+.05),double(yl(2)-1),['P<' num2str(p_high_w)],'Fontsize',14)
set(gca,'TickDir','out','FontSize',14), box off
if isempty(dims_th)
    xlabel(['Weights all modes'])
else
    xlabel(['Weights modes 1:' num2str(dims_th)])
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
             if length(dims) > 1
                w1 = datasets{d}.dim_red_FR{t1}.w(:,1:dims_scatter);
                w2 = datasets{d}.dim_red_FR{t1}.w(:,1:dims_scatter);
             elseif dims < 0
                w1 = datasets{d}.dim_red_FR{t1}.w(:,end+dims_scatter+1:end);
                w2 = datasets{d}.dim_red_FR{t2}.w(:,end+dims_scatter+1:end);
             end
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


% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% % Scatter plot of all pairs of weights for all modes in all task
% % comparisons
% % --Note that this figure compares the weight onto the i-th mode in task 1
% % to the weight onto the i-th mode in task 2. So it looks cool but it isn't
% % very informative
% 
% col_wrist = [247 152 80]/255;
% col_reach = [157 157 182]/255;
% 
% % This one doesn't take the absolute value
% figure,hold on
% plot(all_w_pairs(all_w_pairs_wrist_flg,1),all_w_pairs(all_w_pairs_wrist_flg,2),...
%     '.','linestyle','none','color',col_wrist)
% plot(all_w_pairs(~all_w_pairs_wrist_flg,1),all_w_pairs(~all_w_pairs_wrist_flg,2),...
%     '.','linestyle','none','color',col_reach)
% % high-weight threshold
% plot([-1 -w_th],[w_th w_th],'r','linewidth',2)
% plot([w_th 1],[w_th w_th],'r','linewidth',2)
% text(-.9,.9,'High weight','color','r','FontSize',12)
% plot([-1 -w_th],[-w_th -w_th],'r','linewidth',2)
% plot([w_th 1],[-w_th -w_th],'r','linewidth',2)
% plot([-w_th -w_th],[-1 -w_th],'r','linewidth',2)
% plot([-w_th -w_th],[w_th 1],'r','linewidth',2)
% plot([w_th w_th],[-1 -w_th],'r','linewidth',2)
% plot([w_th w_th],[w_th 1],'r','linewidth',2)
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% xlim([-1 1]),ylim([-1 1])
% xlabel('Weight task 1'),ylabel('Weight task 2')
% if isempty(dims_scatter)
%     title('PCA weights across tasks for all neural modes')
% else
%     title(['PCA weights across tasks for all the leading ' num2str(dims_scatter(end)) ' modes'])
% end
% 
% 
% % And this one does take the absolute value
% figure,hold on
% plot(abs(all_w_pairs(all_w_pairs_wrist_flg,1)),abs(all_w_pairs(all_w_pairs_wrist_flg,2)),...
%     '.','linestyle','none','color',col_wrist)
% plot(abs(all_w_pairs(~all_w_pairs_wrist_flg,1)),abs(all_w_pairs(~all_w_pairs_wrist_flg,2)),...
%     '.','linestyle','none','color',col_reach)
% % high-weight threshold
% plot([w_th w_th],[w_th 1],'r','linewidth',2)
% plot([w_th 1],[w_th w_th],'r','linewidth',2)
% text(.8,double(w_th+.05),'High weight','color','r','FontSize',12)
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% xlim([0 1]),ylim([0 1])
% xlabel('Weight task 1'),ylabel('Weight task 2')
% if isempty(dims_scatter)
%     title('PCA weights across tasks for all neural modes')
% else
%     title(['PCA weights across tasks for all the leading ' num2str(dims_scatter(end)) ' modes'])
% end




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
            if length(dims) > 1
                w1 = datasets{d}.dim_red_FR{t1}.w(:,dims);
                w2 = datasets{d}.dim_red_FR{t1}.w(:,dims);
            elseif dims < 0
                w1 = datasets{d}.dim_red_FR{t1}.w(:,end+dims+1:end);
                w2 = datasets{d}.dim_red_FR{t2}.w(:,end+dims+1:end);
            end
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

        % ---------------------
        % high-weight units in both tasks
        if length(uniq_ro)<length(uniq_ro2)
            u_high_w_both = uniq_ro(ismember(uniq_ro,uniq_ro2));
        else
            u_high_w_both = uniq_ro2(ismember(uniq_ro2,uniq_ro));
        end

        % ---------------------
        % units with high weight in one of the tasks
        high_w_units_in1task_only = union(ro,ro2);


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
        all_high_w(ctr).high_w_units_in1task_only = high_w_units_in1task_only;
        all_high_w(ctr).n_high_w_units_in1task_only = length(high_w_units_in1task_only);
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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Percentage of units that have a high weight in one task or across two
% tasks, for wrist and reach separately


% Do for both datasets combined
n_high_w_units_all_ds = sum([all_high_w.n_high_w_units]);
n_units_all_ds = sum([all_high_w.n_units]); % Note that this is twice the number because we compare 1 to 2 and 2 to 1
n_high_w_units_only1task_all_ds = sum([all_high_w.n_high_w_units_in1task_only]);

perc_high_w_both = n_high_w_units_all_ds/n_high_w_units_only1task_all_ds*100;

% Do separately: for wrist
n_high_w_units_all_wrist = sum([all_high_w([all_high_w.wrist_ds]).n_high_w_units]);
n_units_all_wrist = sum([all_high_w([all_high_w.wrist_ds]).n_units]); % Note that this is twice the number because we compare 1 to 2 and 2 to 1
n_high_w_units_only1task_all_wrist = sum([all_high_w([all_high_w.wrist_ds]).n_high_w_units_in1task_only]);

perc_high_w_both_wrist = n_high_w_units_all_wrist/n_high_w_units_only1task_all_ds*100;

abs_perc_high_w_both_wrist = n_high_w_units_all_wrist/n_units_all_wrist*100;
rel_perc_high_w_only1_wrist = (n_high_w_units_only1task_all_wrist-n_high_w_units_all_wrist)/n_units_all_wrist*100;


% and for reaching
n_high_w_units_all_reach = sum([all_high_w(~[all_high_w.wrist_ds]).n_high_w_units]);
n_units_all_reach = sum([all_high_w(~[all_high_w.wrist_ds]).n_units]); % Note that this is twice the number because we compare 1 to 2 and 2 to 1
n_high_w_units_only1task_all_reach = sum([all_high_w(~[all_high_w.wrist_ds]).n_high_w_units_in1task_only]);

perc_high_w_both_reach = n_high_w_units_all_reach/n_high_w_units_only1task_all_ds*100;

abs_perc_high_w_both_reach = n_high_w_units_all_reach/n_units_all_reach*100;
rel_perc_high_w_only1_reach = (n_high_w_units_only1task_all_reach-n_high_w_units_all_reach)/n_units_all_reach*100;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% BAR Plots summarizing percentages
% figure, hold on
% bar([perc_high_w_both_wrist 100-perc_high_w_both_wrist; perc_high_w_both_reach 100-perc_high_w_both_reach],'Stacked')
% text(.8,perc_high_w_both_wrist-5,[num2str(perc_high_w_both_wrist,3) ' %'],'Fontsize',12,'color','w')
% text(1.8,perc_high_w_both_reach-5,[num2str(perc_high_w_both_reach,3) ' %'],'Fontsize',12,'color','w')
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% set(gca,'XTick',1:2,'XTickLabel',{'Wrist','Reach-to-grasp'},'XTickLabelRotation',45)
% xlim([0 6])
% legend('Units w high weight 2 tasks','All units'),legend boxoff
% ylabel('Percentage units')


figure, hold on
bar([abs_perc_high_w_both_wrist rel_perc_high_w_only1_wrist 100-(abs_perc_high_w_both_wrist+rel_perc_high_w_only1_wrist); ...
    abs_perc_high_w_both_reach rel_perc_high_w_only1_reach 100-(abs_perc_high_w_both_reach+rel_perc_high_w_only1_reach)],'Stacked')
text(.8,abs_perc_high_w_both_wrist-5,[num2str(abs_perc_high_w_both_wrist,3) ' %'],'Fontsize',12,'color','w')
text(1.8,abs_perc_high_w_both_reach-5,[num2str(abs_perc_high_w_both_reach,3) ' %'],'Fontsize',12,'color','w')
text(.8,abs_perc_high_w_both_wrist+rel_perc_high_w_only1_wrist-5,[num2str(rel_perc_high_w_only1_wrist,3) ' %'],'Fontsize',12,'color','w')
text(1.8,abs_perc_high_w_both_reach+rel_perc_high_w_only1_reach-5,[num2str(rel_perc_high_w_only1_reach,3) ' %'],'Fontsize',12,'color','w')
text(.8,95,[num2str(100-(abs_perc_high_w_both_wrist+rel_perc_high_w_only1_wrist),3) ' %'],'Fontsize',12,'color','w')
text(1.8,95,[num2str(100-(abs_perc_high_w_both_reach+rel_perc_high_w_only1_reach),3) ' %'],'Fontsize',12,'color','w')
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
set(gca,'XTick',1:2,'XTickLabel',{'Wrist','Reach-to-grasp'},'XTickLabelRotation',45)
xlim([0 6])
legend('Units w high weight 2 tasks','Units w high weight 1 task','All units'),legend boxoff
ylabel('Percentage units')





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To address Reviewer 1 Comment 7: Look at whether the leading modes
% reflect only a few neurons or are a population wide phenomenon
%
% The way we did it was compute the histograms of weights onto each mode,
% and then compute the ratio between the number of high weight units and...

% units with high weight to the
% units with low weight. If the ratio were 


x_ax = 0:0.01:1.0;

% thresholds for low weight and high-weight units
th_low_weight_hist = 0.2;
th_high_weight_hist = w_th;


% preallocate
n_low_weight_units = [];
n_high_weight_units = [];
R2_high_weight_unit_mode = [];
R2_low_weight_unit_mode = [];

for d = 1:length(datasets)
    
    tasks = length(datasets{d}.labels);
    
    for t = 1:length(tasks)

        % get weights
        W = abs(datasets{d}.dim_red_FR{t}.w(:,dims));

        % for each mode, get the distribution of weights and the number of
        % high and low weights
        for m = 1:size(W,2)
            
            % distribution of weights
            histW(m,:) = histcounts(W(:,m),x_ax);
            
            % ding number of high and low weight units
            n_low_weight_units = [n_low_weight_units; sum(W(:,m)<th_low_weight_hist)];
            n_high_weight_units = [n_high_weight_units; sum(W(:,m)>th_high_weight_hist)];
           
            % if there are any high-weight units see how correlated they
            % are with the mode dynamics
            if n_high_weight_units(end) > 0
                
                % compute correlation between the activity of that unit and
                % the neural mode
                t_high_weight_u = find(W(:,m)>th_high_weight_hist);
                
                for hw = 1:length(t_high_weight_u)
                    psth = datasets{d}.stdata{t}.target{end}.neural_data.conc_smoothed_fr(:,t_high_weight_u(hw));
                    modedyn = datasets{d}.stdata{t}.target{end}.neural_data.dim_red.scores(:,m);
                    
                    R2_high_weight_unit_mode = [R2_high_weight_unit_mode ; abs(calc_r(detrend(psth),detrend(modedyn)))];
                end
            end
            
            % do for low weight units
            if n_low_weight_units(end) > 0
                
                % compute correlation between the activity of that unit and
                % the neural mode
                t_low_weight_u = find(W(:,m)<=th_low_weight_hist);
                
                for hw = 1:length(t_low_weight_u)
                    psth = datasets{d}.stdata{t}.target{end}.neural_data.conc_smoothed_fr(:,t_low_weight_u(hw));
                    modedyn = datasets{d}.stdata{t}.target{end}.neural_data.dim_red.scores(:,m);
                    
                    R2_low_weight_unit_mode = [R2_low_weight_unit_mode ; abs(calc_r(detrend(psth),detrend(modedyn)))];
                end
            end
        end        
    end
end


% Plot histogram
hist_R2_low_weight = histcounts(R2_low_weight_unit_mode,x_ax);
hist_R2_high_weight = histcounts(R2_high_weight_unit_mode,x_ax);

figure, hold on
bar(x_ax(1:end-1),hist_R2_low_weight,'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
bar(x_ax(1:end-1),hist_R2_high_weight,'FaceColor','r','EdgeColor','r')
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
set(gca,'Yscale','log')
xlim([0 1]),ylabel('Counts'),xlabel('Correlation with neural mode dynamics')
legend(['Low weight (<' num2str(th_low_weight_hist,3) ')'],['High weight (>' num2str(th_high_weight_hist,3) ')'],...
    'Location','NorthEast'), legend boxoff






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Take all the results and put them into a structure in which each field is
% one dataset


% Compute the actual number of high-weight units
% --need to do this because they can be repeated for the different task
% comparisons
for d = 1:length(datasets)
   
    % get elements of 'all_high_w' that correspong to this ds
    f_this_d = find([all_high_w.ds]==d);
    
    % get unique neurons for this ds
    high_u_ds = cell2mat({all_high_w(f_this_d).high_w_units}');
    uniq_high_w_t_ds = unique(high_u_ds);
    n_uniq_high_w_t_ds = length(uniq_high_w_t_ds);
    
        
    % Get the weights for each high-weight unit across all tasks in this
    % dataset
    all_weights_high_w = [];
%    all_weights_high_w_norm_axis = [];
    sum_weights_high_w = zeros(n_uniq_high_w_t_ds,length(dims));
    mn_weights_high_w = zeros(n_uniq_high_w_t_ds,length(dims));
    for u = 1:n_uniq_high_w_t_ds
        
        all_weights_high_w{u} = [];
        
        for j = 1:length(f_this_d)
            
            idx_u = find( all_high_w(f_this_d(j)).high_w_units == uniq_high_w_t_ds(u) );
            t_w = all_high_w(f_this_d(j)).modes1(idx_u,:);
            all_weights_high_w{u} = [all_weights_high_w{u}; t_w ]; 
%            all_weights_high_w_norm_axis{u} = [all_weights_high_w_norm_axis; ];
        end
        sum_weights_high_w(u,:) = sum(abs(all_weights_high_w{u}),1);
        mn_weights_high_w(u,:) = mean(abs(all_weights_high_w{u}),1);
    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    
    % Store
%     high_w_per_ds(d).high_u_ds = high_u_ds;
%     high_w_per_ds(d).n_high_u_ds = length(high_u_ds);
%     high_w_per_ds(d).modes_high_u_ds 
    high_w_per_ds(d).uniq_high_w_t_ds = uniq_high_w_t_ds;
    high_w_per_ds(d).all_weights_high_w = all_weights_high_w;
%     high_w_per_ds(d).all_weights_high_w_norm_axis 
    high_w_per_ds(d).sum_weights_high_w = sum_weights_high_w;
    high_w_per_ds(d).mn_weights_high_w = mn_weights_high_w;
    high_w_per_ds(d).n_uniq_high_w_t_ds = n_uniq_high_w_t_ds;
    high_w_per_ds(d).n_u_total = length(datasets{d}.neural_chs);
    if ismember(d,wrist_ds)
        high_w_per_ds(d).wrist_ds = 1;
    else
        high_w_per_ds(d).wrist_ds = 0;
    end
end



% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % PLOT THE WEIGHTS OF THE HIGH-WEIGHT UNITS ONTO DIFFERENT MODES
% 
% 
% % -------------------------------------------------------------------------
% % plot the weights onto all the modes of the neurons that have >= 1 high
% % weight 
% 
% % dataset to plot
% dspl = 7;
% 
% 
% idx_r = 1;
% 
% figure, hold on
% for u = 1:high_w_per_ds(dspl).n_uniq_high_w_t_ds
%     wpl = high_w_per_ds(dspl).all_weights_high_w{u}; 
%     imagesc('XData',1:high_w_per_ds(dspl).n_u_total,'YData',[idx_r idx_r+size(wpl,1)-1],'CData',abs(wpl))
%     idx_r = idx_r + size(wpl,1);
% end
% t_y = 0.5;
% for u = 1:high_w_per_ds(dspl).n_uniq_high_w_t_ds
%     n_rows = size(high_w_per_ds(dspl).all_weights_high_w{u},1);
%     t_y = t_y + n_rows;
%     plot([0.5 high_w_per_ds(dspl).n_u_total], [t_y t_y],'r','linewidth',2)
% end
% set(gca,'xlim',[0.5 high_w_per_ds(dspl).n_u_total],'ylim',[.5 t_y]);
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% ylabel('Abs(weight) --units separated by red lines')
% xlabel('Neural mode')
% 
% 
% % -------------------------------------------------------------------------
% % plot the mean and sum of the weights of the neurons that have >= 1 high
% % weight
% 
% 
% figure, 
% subplot(121), imagesc(high_w_per_ds(dspl).mn_weights_high_w), colorbar
% xlabel('Mean abs(weights)'), ylabel('High weight unit*')
% set(gca,'TickDir','out','FontSize',14), box off
% subplot(122), imagesc(high_w_per_ds(dspl).sum_weights_high_w), colorbar
% xlabel('Sum abs(weights)'), ylabel('High weight unit*')
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% 
% 
% % -------------------------------------------------------------------------
% % Plot all the weights for each session
% 
% figure, hold on
% cols_this = parula(length(wrist_ds)+1);
% for d = 1:length(wrist_ds)
%     plot(abs(cell2mat(high_w_per_ds(d).all_weights_high_w')'),'color',cols_this(d,:))
% end
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% xlabel('Neural mode'), ylabel('Weight')
% title('Wrist datasets')
% 
% 
% figure, hold on
% cols_this = parula(length(reach_ds)+1);
% for d = 1:length(reach_ds)
%     plot(abs(cell2mat(high_w_per_ds(d).all_weights_high_w')'),'color',cols_this(d,:))
% end
% set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
% xlabel('Neural mode'), ylabel('Weight')
% title('Reach-to-grasp datasets')





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT THE DISTRIBUTION OF WEIGHTS OF THE HIGH NEURON WEIGHTS PER SESSION TYPE


col_wrist = [247 152 80]/255;
col_reach = [157 157 182]/255;



% Compute distribution of high-weights for high-weight units

x_ax = 0:step_hist:1+step_hist;


w_high_w_wrist = [];
for d = 1:length(wrist_ds)
    w_high_w_wrist = [w_high_w_wrist; abs(reshape(cell2mat(high_w_per_ds(wrist_ds(d)).all_weights_high_w'),[],1))];
end
w_high_w_reach = [];
for d = 1:length(reach_ds)
    w_high_w_reach = [w_high_w_reach; abs(reshape(cell2mat(high_w_per_ds(reach_ds(d)).all_weights_high_w'),[],1))];   
end
    
hist_w_high_w_wrist = histcounts(w_high_w_wrist,x_ax);
hist_w_high_w_reach = histcounts(w_high_w_reach,x_ax);


% see percentage high weights for the high-weight units
perc_w_high_w_wrist = sum(hist_w_high_w_wrist(x_ax(1:end-1)>w_th))/sum(hist_w_high_w_wrist)*100;
perc_w_high_w_reach = sum(hist_w_high_w_reach(x_ax(1:end-1)>w_th))/sum(hist_w_high_w_reach)*100;
perc_high_w_all = sum(hist_all_w(x_axis_all_w(1:end-1)>w_th))/sum(hist_all_w)*100;


% Plot

figure
subplot(121),hold on
bar(x_ax(1:end-1),hist_w_high_w_wrist/sum(hist_w_high_w_wrist)*100,...
        'FaceColor',col_wrist,'EdgeColor',col_wrist)
ah = bar(x_axis_all_w(1:end-1),100*hist_all_w/numel(all_w),...
        'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6]);
alpha(ah,0.5);    
set(gca,'TickDir','out','FontSize',14), box off, title('Wrist datasets')
set(gca, 'YScale', 'log'), 
yl = ylim; plot([w_th w_th],yl,'r','linewidth',1.5)
legend('high-weight units','all units','high-weight thres.'), legend boxoff
ylabel('Percentage (%)')
text(0.6,2,['Prop. High weights: ' num2str(perc_w_high_w_wrist/perc_high_w_all,3)],'Fontsize',12)
if isempty(dims)
    xlabel('Weights all modes')
else
    xlabel(['Weights modes 1:' num2str(dims)])
end


subplot(122),hold on
bar(x_ax(1:end-1),hist_w_high_w_reach/sum(hist_w_high_w_reach)*100,'FaceColor',col_reach,'EdgeColor',col_reach)
ah = bar(x_axis_all_w(1:end-1),100*hist_all_w/numel(all_w),...
        'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6]);
alpha(ah,0.5);   
set(gca,'TickDir','out','FontSize',14), box off, title('Wrist datasets')
set(gca, 'YScale', 'log'), 
yl = ylim; plot([w_th w_th],yl,'r','linewidth',1.5)
legend('high-weight units','all units','high-weight thres.'), legend boxoff
ylabel('Percentage (%)')
text(0.6,2,['Prop. High weights: ' num2str(perc_w_high_w_reach/perc_high_w_all,3)],'Fontsize',12)
if isempty(dims)
    xlabel('Weights all modes')
else
    xlabel(['Weights modes 1:' num2str(dims)])
end
set(gcf,'color','w')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % FIND POSITION (MODE TO WHICH THEY CONTRIBUTE) OF THE HIGH-WEIGHT IN THE
% % HIGH-WEIGHT UNITS
% 
% 
% pos_high_w_high_w_unit_wrist = [];
% pos_high_w_high_w_unit_reach = [];
% norm_pos_high_w_high_w_unit_wrist = [];
% norm_pos_high_w_high_w_unit_reach = [];
% 
% for d = 1:length(wrist_ds)
%     % get all weights for this session
%     all_t = cell2mat(high_w_per_ds(wrist_ds(d)).all_weights_high_w');
%     
%     % find position high weights
%     [r, c] = find(abs(all_t)>w_th);
%     
%     
%     % now get the normalized position
% end
% 
% for d = 1:length(reach_ds)
%     
% end