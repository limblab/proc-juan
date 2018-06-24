%
% Use TME to establish the smallness of the principal angles
%


rng('shuffle', 'twister') % randomize the seed


% TME parameters
surrogate_type = 'surrogate-TC';
n_surrogates = 10000;

% Datasets to use?
ds_to_use = 1:11;

% plot for each dataset?
indiv_plot_flg = false;


% define which datasets are wrist and which ones are reach-to-grasp separately
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];


% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% Load our shuffled principal angles
load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets_original_submission.mat')


% Load CCA params ---to use the same analysis window as in the paper
proj_params = batch_compare_manifold_projs_defaults();
% overwrite manifold dimension, if necessary
proj_params.dim_manifold = 12;



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% DO!


% define matrices to store variables
all_actual_PAs = [];
TME_th = [];
our_shuffle_th = [];
wrist_flg = [];
session_nbr = [];
ctr = 1;


% Do !!!
for d = 1:length(ds_to_use)

    
    % idx dataset
    t_ds = ds_to_use(d);
    
    % task combinations
    comb_tasks = nchoosek(1:length(datasets{t_ds}.labels),2);
    
    
    % ---------------------------------------------------------------------
    % PREPARE THE DATA WINDOWS -- Taken from canon_corr_all_manifolds. Has
    % to be done for all the tasks in the session 
    
    stdata = datasets{t_ds}.stdata;
    
    % 1) equalize trial duration across all tasks
    stdata = equalize_single_trial_dur( stdata, ...
        'time_win', proj_params.time_win(t_ds,:) );
    
    % 2) equalize number of trials for all targets of a given task
    for i = 1:length(stdata)
        stdata{i} = equalize_nbr_trials_p_target( stdata{i} );
    end
    
    % 3) equalize number of trials across tasks
    stdata = equalize_nbr_trials_across_tasks( stdata, 'all_conc' );

    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % do for all combinations
    for c = 1:size(comb_tasks,1)
    
        % idx tasks
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);

        
        % PUT THE DATA INTO THE RIGHT FORMAT: time x unit x condition --I think

        n_bins = length(stdata{t1}.target{1}.t);
        n_targets = length(stdata{t1}.target)-1;
        n_units = length(stdata{t1}.target{1}.neural_data.neural_chs);

        fr1 = zeros(n_bins,n_units,n_targets);
        fr2 = zeros(n_bins,n_units,n_targets);

        for g = 1:length(stdata{t1}.target)-1
%             fr1(:,:,g) = stdata{t1}.target{g}.neural_data.dim_red.st_scores_mn;
            fr1(:,:,g) = stdata{t1}.target{g}.neural_data.smoothed_fr_mn;
        end
        for g = 1:length(stdata{t2}.target)-1
%             fr2(:,:,g) = stdata{t2}.target{g}.neural_data.dim_red.st_scores_mn;
            fr2(:,:,g) = stdata{t2}.target{g}.neural_data.smoothed_fr_mn;
        end


        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        %% GET PRIMARY FEATURES OF ORIGINAL DATA


        % Quantify primary features of the original data
        [targetSigmaT1, targetSigmaN1, targetSigmaC1, M1] = extractFeatures_J(fr1);
        [targetSigmaT2, targetSigmaN2, targetSigmaC2, M2] = extractFeatures_J(fr2);


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
        else
            error('Need to code these TME parameters')
        end


        % Fit the maximum entropy distribution
        maxEntropy1 = fitMaxEntropy(params1); 
        maxEntropy2 = fitMaxEntropy(params2); 


        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        %% COMPUTE PRINCIPAL ANGLES BETWEEN SURROGATE AND ACTUAL DATASETS

        surrPAs = zeros(n_surrogates,proj_params.dim_manifold);

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

            [~, w1] = pca(surr_tensor1);
            [~, w2] = pca(surr_tensor2);

            surrPAs(s,:) = principal_angles(w1(:,1:proj_params.dim_manifold),w2(:,1:proj_params.dim_manifold));
        end


        % Compute principal angles between the real datasets
        actualPAs = principal_angles(datasets{t_ds}.dim_red_FR{t1}.w(:,1:proj_params.dim_manifold),...
                                        datasets{t_ds}.dim_red_FR{t2}.w(:,1:proj_params.dim_manifold));

        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        %% STORE DATA AND PLOT

        % Our shuffled control
        our_shuffled_PA = angle_non_orth(:,1,space_dim==n_units)';
        % Get TME surrogate with same significance
        TME_PA = prctile(surrPAs,P_orth*100); 

        
        % define matrices to store variables
        all_actual_PAs = [all_actual_PAs; actualPAs']; %#ok<*AGROW>
        all_TME_PAs{ctr} = surrPAs; %#ok<*SAGROW>
        TME_th = [TME_th; TME_PA];
        our_shuffle_th = [our_shuffle_th; our_shuffled_PA];
        if ismember(t_ds,wrist_ds), wyn = 1; else wyn = 0; end
        wrist_flg = [wrist_flg; wyn];
        session_nbr = [session_nbr; t_ds];

        ctr = ctr + 1;

        if indiv_plot_flg

            % Plot the actual PAs, our shuffled PA and the surrogate PAs
            hf = figure; hold on
            plot(rad2deg(actualPAs'),'b','linewidth',1.5)
            plot(our_shuffled_PA','--r','linewidth',1.5)
            plot(rad2deg(TME_PA'),'--k','linewidth',1.5)
            plot(rad2deg(surrPAs'),'color',[.6 .6 .6])
            legend('actual','our shuffled','TME','Location','SouthEast'), legend boxoff
            set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(hf, 'color', [1 1 1]);
            xlim([0 proj_params.dim_manifold]), ylim([0 90])
            xlabel('Neural mode'),ylabel('Principal angle (deg)')
            pause; close all
        end
    end
end


% -------------------------------------------------------------------------
%% COMPARISON OF OUR SHUFFLING CONTROL AND THE TME CONTROL



% Compare TME and our shuffled threshold Scatter plot if mani_dim is 12

if proj_params.dim_manifold == 12

    lfit = polyfit(rad2deg(TME_th),our_shuffle_th,1);
    xfit4plot = [rad2deg(min(min(min(TME_th),min(our_shuffle_th))))-5 90];
    yfit4plot = polyval(lfit,xfit4plot);

    % compute correlation
    [r, Pr] = corr(reshape(rad2deg(TME_th),[],1),reshape(our_shuffle_th,[],1));

    hf = figure; hold on
    plot([0 90],[0 90],'color',[.6 .6 .6])
    plot(xfit4plot,yfit4plot,'k','linewidth',1.5)
    plot(rad2deg(TME_th),our_shuffle_th,'.k','markersize',12)
    legend('perfectly equal','method match','Location','SouthEast'),legend boxoff
    set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(hf, 'color', [1 1 1]);
    text(10,75,[num2str(lfit(2),3) '+' num2str(lfit(1),3) '·x'],'FontSize',14)
    text(10,65,['r=' num2str(r,3) '; P=' num2str(Pr,3)],'FontSize',14)
    xlabel('TME threshold'); ylabel('Our random sampling threshold')
    xlim([0 90]),ylim([0 90])
end



% -------------------------------------------------------------------------
%%

% Plot histogram number similar modes
diff_TME_actual = TME_th - all_actual_PAs;

small_PAs = diff_TME_actual > 0;

highest_similar_mode = zeros(1,length(session_nbr));

for c = 1:length(session_nbr)
    hsm = find(small_PAs(c,:)==0,1);
    if isempty(hsm)
        highest_similar_mode(c) = size(diff_TME_actual,2);
    else
        highest_similar_mode(c) = hsm-1;
    end
end

[c_sim_modes, x_ax] = histcounts( highest_similar_mode, 1:proj_params.dim_manifold+1 );
c_sim_modes = c_sim_modes / sum(c_sim_modes) * 100;

% % Histogram without specifying the task comparison
% figure,
% hb = bar(x_ax(1:end-1),c_sim_modes,'histc');
% set(hb,'FaceColor',[.6 .6 .6])
% set(gca,'TickDir','out','FontSize',14), box off
% ylabel('Percentage (%)')
% xlabel('Highest similar mode')


% Histogram with all task comparisons, organized per task-comparison type

meta_info = batch_get_monkey_task_data(datasets);


lgn_hist        = cell(1,length(meta_info.task_pairs.unique_pairs));
for p = 1:length(meta_info.task_pairs.unique_pairs)
    lgn_hist{p} = [meta_info.task_pairs.unique_pairs{p}{1} ' vs. ' ...
                    meta_info.task_pairs.unique_pairs{p}{2}];
end


highest_similar_mode_type = zeros(length(meta_info.task_pairs.unique_pairs),proj_params.dim_manifold);

for i = 1:length(meta_info.task_pairs.task_pair_nbr)
    
    thm = highest_similar_mode(i);
    ttp = meta_info.task_pairs.task_pair_nbr(i);
    
    highest_similar_mode_type(ttp,thm) = highest_similar_mode_type(ttp,thm) + 1;
end

% convert to %
highest_similar_mode_type = highest_similar_mode_type./sum(sum(highest_similar_mode_type))*100;
highest_similar_mode_type(isnan(highest_similar_mode_type )) = 0;


% -------------------------------------------------------------------------
%% 	PAPER PLOTS (FIGURE 3)


% plot specific sessions --won't do anything if empty
sessions2plot = [11, 8];


% Plot specific sessions --actual PAs, our shuffled PA and the surrogate PAs 

if ~isempty(sessions2plot)
    
    for p = 1:length(sessions2plot)
       
        idx_this = find(session_nbr==sessions2plot(p));
        colors = parula(length(idx_this));
        if size(colors,1) == 1, colors = [0 0 0]; end
        
        hf = figure; hold on
        plot(our_shuffle_th(idx_this(1),:)','--','color',[.5 .5 .5],'linewidth',1.5)
        plot(rad2deg(TME_th(idx_this(1),:)'),'--','color',[.5 0 .5],'linewidth',1.5)
        for c = 1:size(colors,1)
            plot(rad2deg(all_actual_PAs(idx_this(c),:)'),'color',colors(c,:),'linewidth',1.5)
        end
        plot(rad2deg(all_TME_PAs{find(session_nbr==sessions2plot(p),1)}'),'color',[216 190 216]/255)
        legend('our shuffled','TME','Location','SouthEast'), legend boxoff
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(hf, 'color', [1 1 1]);
        xlim([0 proj_params.dim_manifold]), ylim([0 90])
        xlabel('Neural mode'),ylabel('Principal angle (deg)')
        hf.Renderer = 'Painters';
        set(gcf,'color','w')
    end
end





% -------------------------------------------------------------------------
%% PAPER FIGURE 3: EXAMPLE SESSIONS; SUMMARY HISTOGRAM

figure
% Example sessions
for p = 1:length(sessions2plot)
    
    idx_this = find(session_nbr==sessions2plot(p));
    colors = parula(length(idx_this));
    if size(colors,1) == 1, colors = [0 0 0]; end
    
    subplot(1,length(sessions2plot)+1,p), hold on

    % get most stringent significance threshold
    idx_all_comps_this = find(session_nbr==sessions2plot(p));
    TME_th_this = TME_th(idx_all_comps_this,:);
    if numel(idx_all_comps_this ) > 1
        stringent_TME = min(TME_th_this);
    else
        stringent_TME = TME_th_this;
    end
    
    % plot
    for c = 1:size(colors,1)
        plot(rad2deg(all_actual_PAs(idx_this(c),:)'),'color',colors(c,:),'linewidth',1.5)
    end
    plot(rad2deg(stringent_TME'),'--','color',[.6 .6 .6],'linewidth',1.5)
    
    set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off
    xlim([0 proj_params.dim_manifold]), ylim([0 90])
    xlabel('Neural mode')
    if p == 1, ylabel('Principal angle (deg)'); end
end
% Histogram
subplot(1,length(sessions2plot)+1,length(sessions2plot)+1)
bar(x_ax(1:end-1),highest_similar_mode_type','stack')
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off
xlabel('Highest similar mode'),ylabel('Percentage (%)')
legend(lgn_hist,'Location','NorthWest'), legend boxoff
xlim([0 proj_params.dim_manifold+1])
set(gcf,'color','w')


% -------------------------------------------------------------------------
%% SUPPL FIGURE: NORMALIZED PAs

norm_PAs = all_actual_PAs./TME_th;

cols_norm = parula(length(meta_info.task_pairs.unique_pairs));

lgn_ctr = 1;

figure,hold on
for c = 1:size(norm_PAs,1)
    
    ttc = meta_info.task_pairs.task_pair_nbr(c);
    t_c = cols_norm(ttc,:);
    
    if ttc == lgn_ctr
        hp(lgn_ctr)  = plot(norm_PAs(c,:),'color',t_c,'linewidth',1.5);
        lgn_ctr = lgn_ctr + 1;
    else
        plot(norm_PAs(c,:),'color',t_c,'linewidth',1.5);
    end
end
plot([0 proj_params.dim_manifold],[1 1],'--','color',[.6 .6 .6],'linewidth',1.5)
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off
xlabel('Highest similar mode'),ylabel('Normalized principal angle')
ylim([0 1.1]),xlim([0 proj_params.dim_manifold+2])
legend(hp, lgn_hist, 'location','SouthEast'), legend boxoff