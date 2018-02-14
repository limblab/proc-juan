%
% New control for CCA where we shuffle the weights of the linear
% combinations
% 


gcp;

% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% Number of shuffles
n_shuffles = 2000;
P_th = 0.01;

% Plot for each comparison? 
plot_per_comp_flg = true;

% Load params for CCA --- these are the same ones used in the original
% submission of the paper
proj_params = batch_compare_manifold_projs_defaults();

% overwrite manifold dimension, if necessary
proj_params.dim_manifold = 12;

% define which datasets are wrist and which ones are reach-to-grasp,
% separately
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define matrices to store the results
signif_CC.all_actual_CCs = [];
signif_CC.signif_th = [];
signif_CC.wrist_flg = [];
signif_CC.session_nbr = [];
ctr = 1;


% DO

for ds = 1:length(datasets)

    
    % Get task comparisons for this session
    comb_tasks = nchoosek(1:length(datasets{ds}.labels),2);

    
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

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-allocate matrix for saving the results


    % TODO
    
    % Do for each task comparison
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
        
        % Compute the CC
        [A,B,r,U,V] = canoncorr(ssc1,ssc2);


        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SHUFFLE THE WEIGHTS OF THE PCA MODE TO GENERATE DIFFERENT NEUERAL
        % MODES THAT WE'LL COMPARE WITH CCA
        
        % get single-trial neuron activity patterns --matrices are time x
        % neurons x trials 
        fr1 = stdata{comb_tasks(c,1)}.target{end}.neural_data.dim_red.st_scores;
        fr2 = stdata{comb_tasks(c,2)}.target{end}.neural_data.dim_red.st_scores;

        % turn fri into (time x trials) x neurons matrix
        pfr1 = permute(fr1,[1 3 2]);
        sfr1 = reshape(pfr1,[],size(pfr1,3));

        pfr2 = permute(fr2,[1 3 2]);
        sfr2 = reshape(pfr2,[],size(pfr2,3));
        
        
        % variable to store the CCA after projecting onto the shuffled
        % manifolds
        all_rshuff = zeros(n_shuffles,proj_params.dim_manifold);
        
        for s = 1:n_shuffles
        
            % do PCA
            w1 = pca(sfr1);
            w2 = pca(sfr2);
            
            % keep the leading n components
            w1 = w1(:,1:proj_params.dim_manifold);
            w2 = w2(:,1:proj_params.dim_manifold);
            
            % shuffle the weights in the PCA matrices
            idx_shuffle = datasample(1:numel(w1),numel(w1),'Replace',false);
            idx_shuffle = reshape(idx_shuffle,size(w1,1),size(w1,2));
            
            % shuffle the eigenvectors
            w1_shuff = w1(idx_shuffle);
            w2_shuff = w2(idx_shuffle);
            
            % project the neural data onto the manifold
            sc1 = sfr1*w1_shuff;
            sc2 = sfr2*w2_shuff;
            
            
            % Compute the CC
            [~,~,all_rshuff(s,:)] = canoncorr(sc1,sc2);
        end
        
        
        % -----------------------------------------------------------------
        % Store results
        
        signif_CC.all_actual_CCs = [signif_CC.all_actual_CCs; r];
        % Significance threshold
        t_signif_th = prctile(all_rshuff,(1-P_th)*100);
        signif_CC.signif_th = [signif_CC.signif_th; t_signif_th]; %#ok<*AGROW>
        signif_CC.all_shuff_CCs{ctr} = all_rshuff; %#ok<*SAGROW>
        if ismember(ds,wrist_ds), wyn = 1; else wyn = 0; end
        signif_CC.wrist_flg = [signif_CC.wrist_flg; wyn];
        signif_CC.session_nbr = [signif_CC.session_nbr; ds];
        
        ctr = ctr + 1;
        
        % -----------------------------------------------------------------
        % plot
        
        if plot_per_comp_flg
            figure,hold on %#ok<UNRCH>
            plot(r,'k','linewidth',1.5)
            plot(t_signif_th,'color',[.5 0 .5],'linewidth',1.5)
            plot(all_rshuff','color',[216 190 216]/255)
            ylim([0 1]),xlim([0 proj_params.dim_manifold]);
            set(gca,'TickDir','out','FontSize',14), box off
            xlabel('Neural mode'),ylabel('Canonical correlation')
            legend('actual','threshold','surrogates')
            % pause; close
        end
    end
end


% -------------------------------------------------------------------------
%% Summary calculations


% Find number of modes for each task whose CC is above this significance
% threshold

cc_diff = signif_CC.all_actual_CCs-signif_CC.signif_th;

for c = 1:size(cc_diff,1)
    
    highest_similar_mode(c) = find(cc_diff(c,:)<0,1)-1;
end


% Do histogram of the highest modes
x_hist = 1:proj_params.dim_manifold+1;
y_hist = histcounts(highest_similar_mode,x_hist);


% Plot
figure,
bar(x_hist(1:end-1),y_hist,'FaceColor',[.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Highest mode w. similar dynamics'),ylabel('Counts')