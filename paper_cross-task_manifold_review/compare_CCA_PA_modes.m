%
% Compare CCA and PA modes. To answer Reviewr 3 Minor comment 3
%


% load data
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% Load params for CCA --- these are the same ones used in the original
% submission of the paper
proj_params = batch_compare_manifold_projs_defaults();

% overwrite manifold dimension, if necessary
proj_params.dim_manifold = 12;
dims = proj_params.dim_manifold;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare by computing angles between directions


% pre-allocate for saving results
all_dots = []; % dot product between corresponding pairs of CCA and PA modes
all_largest_dots = []; % largest do product between each CCA mode and any PA mode
all_modes_largest_dots = []; % PA mode that yields the previous largest dot product


for ds = 1:length(datasets)
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPARE THE SINGLE TRIALS -- Taken from canon_corr_all_manifolds. Has
    % to be done for all the tasks in the session
    
    stdata = datasets{ds}.stdata;
    
    % 1) equalize trial duration across all tasks
    stdata = equalize_single_trial_dur( stdata, ...
        'time_win', proj_params.time_win(ds,:) );
    
    % 2) equalize number of trials for all targets of a given task
    for i = 1:length(stdata)
        stdata{i} = equalize_nbr_trials_p_target( stdata{i} );
    end
    
    % 3) equalize number of trials across tasks
    stdata = equalize_nbr_trials_across_tasks( stdata, 'all_conc' );

    
    % some definitions
    n_tasks = length(datasets{ds}.labels);
    comb_tasks = nchoosek(1:n_tasks,2);
    
    
    
    for c = 1:size(comb_tasks,1)
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PREPROCESS AND COMPUTE CCs OF THE ACTUAL DATA

        % get single-trial scores -matrices are time x neurons x trials
        sc1 = stdata{comb_tasks(c,1)}.target{end}.neural_data.dim_red.st_scores;
        sc2 = stdata{comb_tasks(c,2)}.target{end}.neural_data.dim_red.st_scores;

        % turn this into (time x trials) x neurons matrix
        psc1 = permute(sc1,[1 3 2]);
        ssc1 = reshape(psc1,[],size(psc1,3));

        psc2 = permute(sc2,[1 3 2]);
        ssc2 = reshape(psc2,[],size(psc2,3));

        % keep only the relevant manifold dimensions
        ssc1 = ssc1(:,1:dims);
        ssc2 = ssc2(:,1:dims);
        
        % Compute the CC matrices A, B
        [A,B,r] = canoncorr(ssc1,ssc2);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE PRINCIPAL ANGLES
        
        W1 = datasets{ds}.dim_red_FR{comb_tasks(c,1)}.w(:,1:dims);
        W2 = datasets{ds}.dim_red_FR{comb_tasks(c,2)}.w(:,1:dims);
        
        % Compute the PA matrices U,V
        [~, U, ~, V]  = principal_angles(W1,W2);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPARE CCA MODES AND PA MODES
        % which is the same as comparing the rows of U and A, V and B
        
        mode1_comp = dot(A/norm(A),U);
        mode2_comp = dot(B/norm(B),V');
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A SMARTER COMPARISON: FIND THE PA MODES THAT HAVE THE MAXIMUM DOT
        % PRODUCT WITH EACH CCA MODE
        
        largest_o1 = zeros(1,dims); % PA mode with largest dot prod with each CCA mode
        largest_dot1 = zeros(1,dims); % largest dot prod of each CCA mode with any PA mode
        largest_o2 = zeros(1,dims); 
        largest_dot2 = zeros(1,dims); 
        
        for mt = 1:dims
            dot_o1 = zeros(1,dims);
            dot_o2 = zeros(1,dims);
            for mo = 1:dims
               dot_o1(mo) = dot(A(:,mt)/norm(A(:,mt)),U(:,mo)); 
               dot_o2(mo) = dot(B(:,mt)/norm(B(:,mt)),V(mo,:)); % because the basis is V' not V
            end
            [largest_dot1(mt), largest_o1(mt)] = max(abs(dot_o1));
            [largest_dot2(mt), largest_o2(mt)] = max(abs(dot_o2));
        end
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STORE VALUES
        
        all_dots = [all_dots; mode1_comp; mode2_comp]; %#ok<*AGROW>
        all_largest_dots = [all_largest_dots; largest_dot1; largest_dot2];
        all_modes_largest_dots = [all_modes_largest_dots; largest_o1; largest_o1];
    end
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS


% -------------------------------
% PLOT pairwise dot products

% figure
% plot(all_dots','color',[.6 .6 .6])


% Compute histograms dot product per dimension

x_hist = 0:0.05:1.05;
hist_dot = zeros(dims,length(x_hist)-1);
for d = 1:dims
    hist_dot(d,:) = histcounts(abs(all_dots(:,d)),x_hist); 
end


% One histogram per dimension --hardcoded for 12
cols_hist = parula(dims+1);

n_rows = 3;
n_cols = 4;

figure('units','normalized','outerposition',[0 0 1 1])
for s = 1:dims
    subplot(n_rows,n_cols,s)
    bar(x_hist(1:end-1),hist_dot(s,:),'FaceColor',cols_hist(s,:))
    set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
    legend(['Mode ' num2str(s)]), legend boxoff
    
    if s > (n_rows-1)*n_cols, xlabel('Cos. angle CCA-PA modes'); end
    if ismember(s,[1 5 9]), ylabel('Counts'); end
end

% histogram pooling all dimensions
mn_all = mean(reshape(abs(all_dots),[],1));
sd_all = std(reshape(abs(all_dots),[],1));

hist_dot_all = sum(hist_dot,1);

figure, hold on
bar(x_hist(1:end-1),hist_dot_all,'FaceColor', [.6 .6 .6], 'EdgeColor', [.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
yl = ylim;
y_stats = (yl(2)-max(hist_dot_all))/2+max(hist_dot_all);
plot(mn_all,y_stats,'.','markersize',32,'color',[.6 .6 .6])
plot([mn_all-sd_all, mn_all+sd_all],[y_stats, y_stats],'color',[.6 .6 .6],'linewidth',1.5)
ylabel('Counts'),xlabel('Cos. angle CCA-PA modes')



% -------------------------------
% PLOT dot largest dot product between each CCA mode and any PA mode


% Compute histograms dot product per dimension

hist_dot_largest = zeros(dims,length(x_hist)-1);
for d = 1:dims
    hist_dot_largest(d,:) = histcounts(abs(all_largest_dots(:,d)),x_hist); 
end


% histogram pooling all dimensions
mn_all_largest = mean(reshape(abs(all_largest_dots),[],1));
sd_all_largest = std(reshape(abs(all_largest_dots),[],1));

hist_dot_all_largest = sum(hist_dot_largest,1);

figure, hold on
bar(x_hist(1:end-1),hist_dot_all_largest,'FaceColor', [.6 .6 .6], 'EdgeColor', [.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
yl = ylim;
y_stats_largest = (yl(2)-max(hist_dot_all_largest))/2+max(hist_dot_all_largest);
plot(mn_all_largest,y_stats_largest,'.','markersize',32,'color',[.6 .6 .6])
plot([mn_all_largest-sd_all_largest, mn_all_largest+sd_all_largest],[y_stats_largest, y_stats_largest],'color',[.6 .6 .6],'linewidth',1.5)
ylabel('Counts'),xlabel('Cos. angle closest CCA-PA modes')

