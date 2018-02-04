%
% Try to understand the relationship between dPCA and PCA modes
% 


% choose manifold dimensionality
manifold_dim = 12;

% Use all wrist datasets for dPCA
dPCA_datasets = [1 2 3 7 8 9];

% What dPCA datasets correspond to each monkey
jaco_datasets = 1:3;
jango_datasets = 4:6;

% Do per-session plots
plots_p_session = true;


% Do dPCA if the data aren't available in the WS
if ~exist('dPCA_results','var')
    
    for i = 1:length(dPCA_datasets)
        dPCA_results{i} = call_dPCA( datasets{dPCA_datasets(i)}.stdata, manifold_dim, false );
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we need to get the neural activity as we have projected it onto
% the manifold -- this code is copied directly from call_dPCA.m


for ds = 1:length(dPCA_datasets) % with respect to dPCA_datasets (!!!)

    
    % data for this session
    stdata = datasets{dPCA_datasets(ds)}.stdata;

    % neural units to use
    neural_chs = stdata{1}.target{1}.neural_data.neural_chs;

    % check if the dimensions in single_trial_data are consistent
    if size(stdata{1}.target{1}.neural_data.smoothed_fr,2) ~= numel(neural_chs)
        error('single trial data does not include all the neural channels')
    end

    % make the target averaged responses for each task have equal length.
    % This is not done in single_trial_analysis.m, where single trial
    % duration is only equalized for each task
    stdata = equalize_single_trial_dur( stdata );

    % get rid of the last target, which is all the concatenated targets
    for i = 1:length(stdata)
        stdata{i}.target(end) = [];
    end


    % ------------------------------------------------------------------------
    % Arrange the data (as described in dpca_demo)

    % N is the number of neurons
    % S is the number of conditions --> tasks in our case
    % D is the number of decisions --> targets in our case
    %       ToDo: so far we are choosing the min, but see if they can be different for each task
    % T is the number of time points --each trial should have the same duration
    % in time !!!
    N                   = numel(neural_chs);
    S                   = numel(stdata);
    D                   = min(cellfun(@(x) length(x.target), stdata));
    T                   = size(stdata{1}.target{1}.neural_data.fr,1);
    % max number of repetitions
    max_trial_num       = 0;
    for i = 1:S
        if max(cellfun(@(x) size(x.neural_data.fr,3), stdata{1}.target )) > max_trial_num
            max_trial_num = max(cellfun(@(x) size(x.neural_data.fr,3), stdata{1}.target ));
        end
    end


    % trial_num: N x S x D
    trial_num           = zeros(N,S,D);


    % firing_rates: N x S x D x T x max_trial_num -- populated with our
    % single_trial_data
    firing_rates        = nan(N,S,D,T,max_trial_num);
    % ·········································································
    % In the 1D tasks, targets 1 to 6 go from left to right, but in the 2D
    % tasks targets are ordered as follows: 5, 7, 8, 6, 4, 2, 1, 3 --beginning
    % at 12 o'clock and going clockwise. They will be paired as 1D/2D: 1/1,
    % 2/2, 3/3, 4/6, 5/7, 6, 8 (as defined in target_order)
    % ·········································································
    target_order        = [1 2 3 4 5 6; 1 2 3 6 7 8]; % row 1: 1D task; row 2: 2D task

    for n = 1:N
        for s = 1:S
            for d = 1:D
                if length(stdata{s}.target) == D
                    tgt = target_order(1,d);
                elseif length(stdata{s}.target) >= D
                    tgt = target_order(2,d);
                end
                trials_this = size(stdata{s}.target{tgt}.neural_data.smoothed_fr(:,n,:),3);
                firing_rates(n,s,d,:,1:trials_this) = ...
                    squeeze(stdata{s}.target{tgt}.neural_data.smoothed_fr(:,n,:));
            end
        end
    end


    % firing_rates_average: N x S x D x T -- these are PSTHs
    firing_rates_avg    = nanmean(firing_rates, 5);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---------------------------------------------------------------------
    % Now do PCA jointly on all the tasks in the session -- I will use the
    % trial averaged data for this
    
    % Prepare the data -- Need a neurons X time matrix
    fra = permute(firing_rates_avg,[1 4 3 2]);
    fra = reshape(fra,size(fra,1),size(fra,2)*size(fra,3)*size(fra,4));
    
    [w_pca, lv_pca] = pca(fra'); % checked that this gives the same as Machens' code
    
    % ---------------------------------------------------------------------
    % 1. Compare PCA and dPCA manifolds using principal angles
    
    PA_pca_dpca = principal_angles( w_pca(:,1:manifold_dim), dPCA_results{ds}.W );
    
    
    % ---------------------------------------------------------------------
    % 2. Do CCA between the PCA and dPCA latent variables --this is for the
    % PCA manifold that spans all tasks
    
    % prepare the dPCA latent variables
    lv_dpca = permute(dPCA_results{ds}.lat_vars_mn,[1 4 3 2]);
    lv_dpca = reshape(lv_dpca,size(lv_dpca,1),size(lv_dpca,2)*size(lv_dpca,3)*size(lv_dpca,4));
    
    [~,~,CC,~,~] = canoncorr(lv_pca(:,1:manifold_dim),lv_dpca');
    

    
    % ---------------------------------------------------------------------
    % 4. Do CCA between the dPCA latent variables and the PCA latent
    % variables for the manifold for each specific task
    % Plus a few other things
    
    % pre-allocate vars
    tasks_this = length(datasets{dPCA_datasets(ds)}.labels);
    CC_single_task = zeros(manifold_dim,tasks_this);
    PA_trial_avg_and_single_task = zeros(manifold_dim,tasks_this);
    PA_single_task_vs_multi_task_PCA = zeros(manifold_dim,tasks_this);
    
    for t = 1:tasks_this
        
        % Prepare the spikes
        fra = permute(firing_rates_avg,[1 4 3 2]);
        fra = squeeze(fra(:,:,:,t));
        fra = reshape(fra,size(fra,1),size(fra,2)*size(fra,3));
        
        % And do PCA for this task
        [w_pca_t, lv_pca_t, eig_t] = pca(fra');

        
        % -----------------------------------------------------------------
        % Compare to the non-trial-averaged manifold
        % -- this is more a sanity check than anything
        PA_tavg_st = principal_angles( datasets{dPCA_datasets(ds)}.dim_red_FR{t}.w(:,1:manifold_dim), ...
                        w_pca_t(:,1:manifold_dim));
          
                    
        % -----------------------------------------------------------------            
        % Compare the multi-task and task-specific PCA manifolds
        PA_spec_multi = principal_angles( w_pca(:,1:manifold_dim), ...
                                            w_pca_t(:,1:manifold_dim));
                    
        
        % -----------------------------------------------------------------
        % Do CCA between dPCA latent variables and task-specific PCA latent
        % variables
        
        % Prepare the dPCA latent variables
        lv_dpca = permute(dPCA_results{ds}.lat_vars_mn,[1 4 3 2]);
        lv_dpca = squeeze(lv_dpca(:,:,:,t));
        lv_dpca = reshape(lv_dpca,size(lv_dpca,1),size(lv_dpca,2)*size(lv_dpca,3));

        % Compare the dPCA and PCA (for this task-specific manifold) latent
        % variables
        [~,~,CC_st,~,~] = canoncorr(lv_pca_t(:,1:manifold_dim),lv_dpca');
        
        % Save the results
        CC_single_task(:,t) = CC_st';
        PA_trial_avg_and_single_task(:,t) = PA_tavg_st;
        PA_single_task_vs_multi_task_PCA(:,t) = PA_spec_multi;
    end
    
    % ---------------------------------------------------------------------
    % Save the results of this comparison

    dPCA_PCA_comp(ds).PA = PA_pca_dpca;
    dPCA_PCA_comp(ds).CC = CC'; 
    dPCA_PCA_comp(ds).CC_single_task = CC_single_task;
    dPCA_PCA_comp(ds).PA_trial_avg_and_single_task = PA_trial_avg_and_single_task;
    dPCA_PCA_comp(ds).PA_task_spec_vs_multitask_PCA = PA_single_task_vs_multi_task_PCA;
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTS



% 1. PAs
figure,
plot(rad2deg([dPCA_PCA_comp.PA]),'linewidth',1.5)
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlim([0 manifold_dim]),ylim([0 90])
ylabel('Principal angles (deg)')
xlabel('Neural mode')
title('Comparison between PCA and dPCA manifolds for all wrist tasks')


% 2. CCs of multi-task PCA manifold and dPCA manifold
figure,hold on
plot([dPCA_PCA_comp.CC],'color',[.6 .6 .6])
plot(mean([dPCA_PCA_comp.CC],2),'color','k','linewidth',1.5)
plot(mean([dPCA_PCA_comp.CC],2)+std([dPCA_PCA_comp.CC],0,2),'-.k','linewidth',1.5)
plot(mean([dPCA_PCA_comp.CC],2)-std([dPCA_PCA_comp.CC],0,2),'-.k','linewidth',1.5)
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlim([0 manifold_dim]),ylim([0 1])
ylabel('CC neural mode dynamics')
xlabel('Neural mode')
title('Comparison between PCA and dPCA manifolds for all wrist tasks')


% 3. Sanity check: PA single-trial and trial-average PCA manifolds
figure
plot(rad2deg([dPCA_PCA_comp.PA_trial_avg_and_single_task]),'color',[.6 .6 .6])
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlim([0 manifold_dim]),ylim([0 90])
ylabel('Principal angles (deg)')
xlabel('Neural mode')
title('Comparison between single-trial and trial-averaged PCA manifolds')


% 4. PA between single-task and multi-task PCA manifolds
figure, hold on
plot(rad2deg([dPCA_PCA_comp.PA_task_spec_vs_multitask_PCA]),'color',[.6 .6 .6])
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlim([0 manifold_dim]),ylim([0 90])
ylabel('Principal angles (deg)')
xlabel('Neural mode')
title('Comparison between task-specific and multi-task PCA manifolds')


% 5. CCs of single-task PCA manifold and dPCA manifold
figure,hold on
plot([dPCA_PCA_comp.CC_single_task],'color',[.6 .6 .6])
plot(mean([dPCA_PCA_comp.CC_single_task],2),'color','k','linewidth',1.5)
plot(mean([dPCA_PCA_comp.CC_single_task],2)+std([dPCA_PCA_comp.CC_single_task],0,2),'-.k','linewidth',1.5)
plot(mean([dPCA_PCA_comp.CC_single_task],2)-std([dPCA_PCA_comp.CC_single_task],0,2),'-.k','linewidth',1.5)
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlim([0 manifold_dim]),ylim([0 1])
ylabel('CC neural mode dynamics')
xlabel('Neural mode')
title('Comparison between task-specific PCA and dPCA manifolds')
