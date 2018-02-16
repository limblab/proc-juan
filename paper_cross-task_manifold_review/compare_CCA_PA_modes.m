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


% pre-allocate for saving results
all_dots = [];

% Compare by computing angles between directions
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
        ssc1 = ssc1(:,1:proj_params.dim_manifold);
        ssc2 = ssc2(:,1:proj_params.dim_manifold);
        
        % Compute the CC matrices A, B
        [A,B,r] = canoncorr(ssc1,ssc2);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE PRINCIPAL ANGLES
        
        W1 = datasets{ds}.dim_red_FR{comb_tasks(c,1)}.w(:,1:proj_params.dim_manifold);
        W2 = datasets{ds}.dim_red_FR{comb_tasks(c,2)}.w(:,1:proj_params.dim_manifold);
        
        % Compute the PA matrices U,V
        [~, U, ~, V]  = principal_angles(W1,W2);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPARE CCA MODES AND PA MODES
        % which is the same as comparing U and A, V and B
        
        mode1_comp = dot(A/norm(A),U);
        mode2_comp = dot(B/norm(B),V');
        
        all_dots = [all_dots; mode1_comp; mode2_comp]; %#ok<*AGROW>
    end
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT dot products

% figure
% plot(all_dots','color',[.6 .6 .6])


% Compute histograms dot product per dimension

x_hist = 0:0.05:1.05;
hist_dot = zeros(proj_params.dim_manifold,length(x_hist)-1);
for d = 1:proj_params.dim_manifold
    hist_dot(d,:) = histcounts(abs(all_dots(:,d)),x_hist); 
end


% One histogram per dimension --hardcoded for 12
cols_hist = parula(proj_params.dim_manifold+1);

n_rows = 3;
n_cols = 4;

figure('units','normalized','outerposition',[0 0 1 1])
for s = 1:proj_params.dim_manifold
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