%
% Calculate Percentage VAF of the neural modes projected with CCA
%


% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% Load params for CCA --- these are the same ones used in the original
% submission of the paper
proj_params = batch_compare_manifold_projs_defaults();

% overwrite manifold dimension, if necessary
proj_params.dim_manifold = 12;

% define which datasets are wrist and which ones are reach-to-grasp,
% separately
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11]; % [4:6 10:11];

% What datasets to use
ds_to_use = [wrist_ds, reach_ds];

% retrieve bin size
bin_size = round(mean(diff(datasets{1}.binned_data{1}.timeframe))*100)/100;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DO

for ds = 1:length(ds_to_use)

    
    % Get task comparisons for this session
    comb_tasks = nchoosek(1:length(datasets{ds_to_use(ds)}.labels),2);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPARE THE SINGLE TRIALS -- Taken from canon_corr_all_manifolds. Has
    % to be done for all the tasks in the session
    
    stdata = datasets{ds_to_use(ds)}.stdata;
    
    % 1) equalize trial duration across all tasks
    stdata = equalize_single_trial_dur( stdata, ...
        'time_win', proj_params.time_win(ds_to_use(ds),:) );
    
    % 2) equalize number of trials for all targets of a given task
    for i = 1:length(stdata)
        stdata{i} = equalize_nbr_trials_p_target( stdata{i} );
    end
    
    % 3) equalize number of trials across tasks
    stdata = equalize_nbr_trials_across_tasks( stdata, 'all_conc' );

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-allocate matrix for saving the results
    
    vaf_CCA{ds}.vaf_CC1 = zeros(size(comb_tasks,1),proj_params.dim_manifold); %#ok<*SAGROW>
    vaf_CCA{ds}.vaf_CC2 = zeros(size(comb_tasks,1),proj_params.dim_manifold);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do for each task comparison
    for c = 1:size(comb_tasks,1)
    
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % PREPROCESS AND COMPUTE CCs

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

        % % I don't think this makes sense, and it wasn't in the original analysis
        % zssc1 = zscore(detrend(ssc1,'constant'),1);
        % zssc2 = zscore(detrend(ssc2,'constant'),1);
        % --- this doesn't seem to change anything, non-surprisingly


        % Compute the CC
        [A,B,r,U,V] = canoncorr(ssc1,ssc2);


        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % VARIANCE ACCOUNTED FOR BY THE CCS
        % 
        % This method deals with the CCA axes not being orthogonal, and is
        % largely based on Machens' method to compute the VAF of each dPC,
        % adapted by Sara
        
        % PCA weights
        w_pca1 = stdata{comb_tasks(c,1)}.target{end}.neural_data.dim_red.w(:,1:proj_params.dim_manifold);
        w_pca2 = stdata{comb_tasks(c,2)}.target{end}.neural_data.dim_red.w(:,1:proj_params.dim_manifold);
        
        
        % get matrices Wi that go from neural space to the CCA axes
        % (passing by the PCA axes that define the manifold)
        W1 = w_pca1*A;
        W2 = w_pca2*B;
        
        
        % now compute the matrices V that define the expansion from the PCA
        % axes back onto neural space
        V1 = W1*inv(A'*A);
        V2 = W2*inv(B'*B);
        
        
        % get firing rates into dPCA format % firing_rates_average: Neurons
        % x Targets x Time
        
        n_units = length(stdata{comb_tasks(c,1)}.target{1}.neural_data.neural_chs);
        n_targets = length(stdata{comb_tasks(c,1)}.target)-1;
        n_bins = length(stdata{comb_tasks(c,1)}.target{1}.t);
        
        fr_avg1 = zeros(n_units,n_targets,n_bins);
        fr_avg2 = zeros(n_units,n_targets,n_bins);
        
        for t = 1:length(stdata{comb_tasks(c,1)}.target)-1

            fr_t1 = stdata{comb_tasks(c,1)}.target{t}.neural_data.smoothed_fr_mn;
            fr_t1 = permute(fr_t1,[2 1]);            
            fr_avg1(:,t,:) = fr_t1;
        end
        % do for task 2 in another loop because the number targets may be
        % different
        for t = 1:length(stdata{comb_tasks(c,2)}.target)-1
            
            fr_t2 = stdata{comb_tasks(c,2)}.target{t}.neural_data.smoothed_fr_mn;
            fr_t2 = permute(fr_t2,[2 1]);
            fr_avg2(:,t,:) = fr_t2;          
        end
        
        % -----------------------------------------------------------------
        % Call Machens' dPCA explained variance function

        vaf_CC1 = dpca_explainedVariance(fr_avg1, W1, V1 );
        vaf_CC2 = dpca_explainedVariance(fr_avg2, W2, V2 );
        
        vaf1 = vaf_CC1.cumulativeDPCA;
        vaf2 = vaf_CC2.cumulativeDPCA;
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store the results
        
        vaf_CCA{ds}.vaf_CC1(c,:) = vaf1;
        vaf_CCA{ds}.vaf_CC2(c,:) = vaf2;    
    end    
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POOL DATA ACROSS SESSIONS TOGETHER


% -------------------------------------------------------------------------
% 1. All wrist and reach-to-grasp tasks together

% make matrices with all task comparisons --there does not seem to be a
% difference between sets of tasks
all_vaf_CC1 = cell2mat( cellfun(@(x) x.vaf_CC1, vaf_CCA, 'UniformOutput',false)' );
all_vaf_CC2 = cell2mat( cellfun(@(x) x.vaf_CC2, vaf_CCA, 'UniformOutput',false)' );

% pool all comparisons together
all_vaf_CC = [all_vaf_CC1; all_vaf_CC2];


% -------------------------------------------------------------------------
% 2. Wrist tasks and Reach-to-grasp tasks separately

% define vars
all_vaf_CC1_wrist = [];
all_vaf_CC2_wrist = [];
all_vaf_CC_wrist = [];

all_vaf_CC1_reach = [];
all_vaf_CC2_reach = [];
all_vaf_CC_reach = [];


for d = 1:length(wrist_ds)
    
    all_vaf_CC1_wrist = [all_vaf_CC1_wrist; vaf_CCA{wrist_ds(d)}.vaf_CC1]; %#ok<*AGROW>
    all_vaf_CC2_wrist = [all_vaf_CC2_wrist; vaf_CCA{wrist_ds(d)}.vaf_CC2];
end


for d = 1:length(reach_ds)
    
    all_vaf_CC1_reach = [all_vaf_CC1_reach; vaf_CCA{reach_ds(d)}.vaf_CC1];
    all_vaf_CC2_reach = [all_vaf_CC2_reach; vaf_CCA{reach_ds(d)}.vaf_CC2];
end


% pool all comparisons together
all_vaf_CC_wrist = [all_vaf_CC1_wrist; all_vaf_CC2_wrist];
all_vaf_CC_reach = [all_vaf_CC1_reach; all_vaf_CC2_reach];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT


% -------------------------------------------------------------------------
% WRIST TASKS AND REACH-TO-GRASP TASKS SEPARATELY


% VAF in neural space
figure, hold on
plot(all_vaf_CC1_wrist(1,:),'color',[.7 .7 .7])
plot(all_vaf_CC1_reach(1,:),'color',[1 .7 .3])
for c = 2:size(all_vaf_CC1_wrist,1)
    plot(all_vaf_CC2_wrist(c,:),'color',[.7 .7 .7])
end
for c = 2:size(all_vaf_CC_reach,1)
    plot(all_vaf_CC_reach(c,:),'color',[1 .7 .3])
end
errorbar(1:proj_params.dim_manifold,mean(all_vaf_CC_wrist,1),std(all_vaf_CC_wrist,0,1),...
    'linestyle','none','linewidth',2,'color','k','marker','.','markersize',32)
errorbar(1:proj_params.dim_manifold,mean(all_vaf_CC_reach,1),std(all_vaf_CC_reach,0,1),...
    'linestyle','none','linewidth',2,'color','r','marker','.','markersize',32)
set(gca,'TickDir','out','FontSize',14), box off
ylabel('Normalized variance accounted for (%)')
xlabel('Neural mode')
legend('wrist','reach-to-grasp','Location','SouthEast'), legend boxoff
    


% -------------------------------------------------------------------------
% WRIST TASKS AND REACH-TO-GRASP TASKS SEPARATELY

