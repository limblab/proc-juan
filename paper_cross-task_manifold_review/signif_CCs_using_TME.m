%
% Using TME to assess the significance of TME
% 


clearvars -except datasets , close all;


% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% Number of shuffles
n_surrogates = 1000;
P_th = 0.01;

% use (concatenated) single trials for CCA
single_trial_CCA_flg = false;

% Surrogate type TME
surrogate_type = 'surrogate-N'; % 'surrogate-N' 'surrogate-T', 'surrogate-TC' implemented

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

% Targets of the iso2D tast to compare (exclude top and bottom ones)
tgts_iso2d = [1:3 6:8];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define matrices to store the results
all_actual_CCs = [];
signif_th = [];
wrist_flg = [];
session_nbr = [];
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
    

        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PREPROCESS AND COMPUTE CCs OF THE ACTUAL DATA

        % get data for CCA ---Slightly convoluted because it compares
        % implements fixes for when the number of targets is not the same
        if single_trial_CCA_flg
            % single-trial scores -matrices are time x neurons x trials
            sc1 = stdata{t1}.target{end}.neural_data.dim_red.st_scores;
            sc2 = stdata{t2}.target{end}.neural_data.dim_red.st_scores;
        else
            % trial-averaged scores -matrices are time x neurons x trials
            sc1 = stdata{t1}.target{end}.neural_data.dim_red.st_scores_mn;
            sc2 = stdata{t2}.target{end}.neural_data.dim_red.st_scores_mn;
            
            % when comparing the reach-to-grasp tasks these numbers won't
            % be the same. Or when comparing Iso 2D to any other movement
            % task
            if size(sc1,1) ~= size(sc2,1)
                % for the reach-to-grasp tasks
                if min([length(stdata{t1}.target),length(stdata{t2}.target)]) == 2
                    % just keep the first target
                    sc1 = stdata{t1}.target{1}.neural_data.dim_red.st_scores_mn;
                    sc2 = stdata{t2}.target{1}.neural_data.dim_red.st_scores_mn;
                % when comparing iso 2d to any other task
                else
                    [~,idx_iso2d] = max([length(stdata{t1}.target),length(stdata{t2}.target)]);
                    if idx_iso2d == 2
                        sc2 = [];
                        % keep targets [1:3 6:8]
                        for tg = 1:length(tgts_iso2d)
                            sc2 = [sc2; stdata{t2}.target{tgts_iso2d(tg)}.neural_data.dim_red.st_scores_mn];
                        end
                    else
                        sc1 = [];
                        % keep targets [1:3 6:8]
                        for tg = 1:length(tgts_iso2d)
                            sc1 = [sc1; stdata{t2}.target{tgts_iso2d(tg)}.neural_data.dim_red.st_scores_mn];
                        end
                    end
                end
            end
        end

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
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        %% GET PRIMARY FEATURES OF ORIGINAL DATA

        
        % PUT THE DATA INTO THE RIGHT FORMAT: time x unit x condition --I think
        n_bins = length(stdata{t1}.target{1}.t);
        n_targets = length(stdata{t1}.target)-1;
        n_units = length(stdata{t1}.target{1}.neural_data.neural_chs);
        

        % get trial-averaged FRs -matrices are time x neurons x trials
        
        % this is a bit weird because it accounts for the cases where we
        % don't compare tasks with the same number of targets
        if length(stdata{t1}.target) ~= length(stdata{t2}.target)
            % when comparing ball to grip
            if length(stdata{t1}.target) == 2 || length(stdata{t2}.target) == 2
                sfr1(:,:,1) = stdata{t1}.target{1}.neural_data.smoothed_fr_mn;
                sfr2(:,:,1) = stdata{t2}.target{1}.neural_data.smoothed_fr_mn;
            % when comparing iso 2D to any of the wrist tasks
            elseif length(stdata{t1}.target) == 9 || length(stdata{t2}.target) == 9
                if length(stdata{t1}.target) == 9
                    for g = 1:length(tgts_iso2d)
                        sfr1(:,:,g) = stdata{t1}.target{tgts_iso2d(g)}.neural_data.smoothed_fr_mn;
                        sfr2(:,:,g) = stdata{t2}.target{g}.neural_data.smoothed_fr_mn;
                    end
                elseif length(stdata{t2}.target) == 9
                    for g = 1:length(tgts_iso2d)
                        sfr1(:,:,g) = stdata{t1}.target{g}.neural_data.smoothed_fr_mn;
                        sfr2(:,:,g) = stdata{t2}.target{tgts_iso2d(g)}.neural_data.smoothed_fr_mn;
                    end
                end
            end
        % when both tasks we compare have the same number of targets
        else
            for g = 1:length(stdata{t1}.target)-1
                sfr1(:,:,g) = stdata{t1}.target{g}.neural_data.smoothed_fr_mn;
                sfr2(:,:,g) = stdata{t2}.target{g}.neural_data.smoothed_fr_mn;
            end
        end
%         for g = 1:length(stdata{t2}.target)-1
%             sfr2(:,:,g) = stdata{t2}.target{g}.neural_data.smoothed_fr_mn;
%         end
        
        % Quantify primary features of the original data
        [targetSigmaT1, targetSigmaN1, targetSigmaC1, M1] = extractFeatures_J(sfr1);
        [targetSigmaT2, targetSigmaN2, targetSigmaC2, M2] = extractFeatures_J(sfr2);


        % define surrogate params
        if strcmp(surrogate_type, 'surrogate-TC')
            params1.margCov{1} = targetSigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = targetSigmaC1;
            params1.meanTensor = M1.TC;

            params2.margCov{1} = targetSigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = targetSigmaC2;
            params2.meanTensor = M2.TC;
        
        elseif strcmp(surrogate_type, 'surrogate-T')
            params1.margCov{1} = targetSigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = [];
            params1.meanTensor = M1.T;
            
            params2.margCov{1} = targetSigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = [];
            params2.meanTensor = M2.T;
        elseif strcmp(surrogate_type, 'surrogate-N')
            params1.margCov{1} = targetSigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = [];
            params1.meanTensor = M1.N;
            
            params2.margCov{1} = targetSigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = [];
            params2.meanTensor = M2.N;
        else
            error('Need to code these TME parameters')
        end


        % Fit the maximum entropy distribution
        maxEntropy1 = fitMaxEntropy(params1); 
        maxEntropy2 = fitMaxEntropy(params2); 


        % -----------------------------------------------------------------
        %
        % Do CCA of the surrogate datasets
        
        % variable to store the CCA after projecting onto the shuffled
        % manifolds
        all_rshuff = zeros(n_surrogates,proj_params.dim_manifold);
        
        for s = 1:n_surrogates
        
            fprintf('surrogate %d from %d \n', s, n_surrogates)

            % generate surrogate datasets
            surr_tensor1 = sampleTME(maxEntropy1);
            surr_tensor2 = sampleTME(maxEntropy2);

            % Do PCA OF THE SURROGATE DATASETS
            surr_tensor1 = permute(surr_tensor1,[2 1 3]);
            surr_tensor2 = permute(surr_tensor2,[2 1 3]);

            surr_tensor1 = reshape(surr_tensor1,size(surr_tensor1,1),[]);
            surr_tensor2 = reshape(surr_tensor2,size(surr_tensor2,1),[]);

            [sc1, w1] = pca(surr_tensor1);
            [sc2, w2] = pca(surr_tensor2);
            
            
            % Compute the CC
            [~,~,all_rshuff(s,:)] = canoncorr(sc1(:,1:proj_params.dim_manifold),...
                sc2(:,1:proj_params.dim_manifold));
        end
        
        
        % -----------------------------------------------------------------
        % Store results
        
        all_actual_CCs = [all_actual_CCs; r];
        % Significance threshold
        t_signif_th = prctile(all_rshuff,(1-P_th)*100);
        signif_th = [signif_th; t_signif_th]; %#ok<*AGROW>
        all_shuff_CCs{ctr} = all_rshuff; %#ok<*SAGROW>
        if ismember(ds,wrist_ds), wyn = 1; else wyn = 0; end
        wrist_flg = [wrist_flg; wyn];
        session_nbr = [session_nbr; ds];
        
        ctr = ctr + 1;
        
        % -----------------------------------------------------------------
        % plot
        
        if plot_per_comp_flg
            figure,hold on %#ok<UNRCH>
            plot(r,'k','linewidth',1.5)
            plot(t_signif_th,'color',[.5 0 .5],'linewidth',1.5)
            plot(all_rshuff','color',[216 190 216]/255)
            plot(r,'k','linewidth',1.5)
            ylim([0 1]),xlim([0 proj_params.dim_manifold]);
            set(gca,'TickDir','out','FontSize',14), box off
            xlabel('Neural mode'),ylabel('Canonical correlation')
            legend('actual','surrogate','threshold','Location','SouthWest')
%            pause; close
        end
        
        clear target* sfr*
    end
end


% -------------------------------------------------------------------------
%% Summary calculations


% Find number of modes for each task whose CC is above this significance
% threshold

cc_diff = all_actual_CCs-signif_th;

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