%
% New control for CCA where we shuffle the weights of the linear
% combinations
% 


% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% Number of shuffles
n_shuffles = 100;

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


    
    % Do for each task comparison
    for c = 1:size(comb_tasks,1)
    
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PREPROCESS AND COMPUTE CCs

        % get single-trial neuron activity patterns --matrices are time x
        % neurons x trials 
        
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
        % SHUFFLE THE WEIGHTS (IN COLUMNS) OF THE LINEAR COMBINATION
        % MATRICES A AND B
        
        % variable to store the shuffles
        all_rshuff = zeros(n_shuffles,proj_params.dim_manifold);
        
        for s = 1:n_shuffles
        
            % generate random idx for shuffling
            rng('shuffle');
            idx_shuffleA = datasample(1:numel(A),numel(A),'Replace',true);
            idx_shuffleA = reshape(idx_shuffleA,size(A,1),size(A,2));
            rng('shuffle');
            idx_shuffleB = datasample(1:numel(B),numel(B),'Replace',true);
            idx_shuffleB = reshape(idx_shuffleB,size(B,1),size(B,2));
            
            % shuffle the linear combination weights using the same indexes
            % for both
            Ashuff = A(idx_shuffleA);
            Bshuff = B(idx_shuffleB);
            
            % project the data
            surr1 = (ssc1-repmat(mean(ssc1),[size(ssc1,1),1]))*Ashuff;
            surr2 = (ssc2-repmat(mean(ssc1),[size(ssc2,1),1]))*Bshuff;
        
            % Compute the CC
            [~,~,all_rshuff(s,:)] = canoncorr(surr1,surr2);
            % all_rshuff(s,:) = abs(calc_r(surr1,surr2));
        end
        
        % plot?
        if plot_per_comp_flg
            figure,hold on
            plot(all_rshuff','color',[.6 .6 .6])
            plot(r,'k','linewidth',1.5)
            ylim([0 1]);
            
            figure
            imagesc(all_rshuff);
        end
    end
end