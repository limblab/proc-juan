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
reach_ds = [4:6 10:11];

% retrieve bin size
bin_size = round(mean(diff(datasets{1}.binned_data{1}.timeframe))*100)/100;



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
    
    vaf_CCA{ds}.vaf_CC1 = zeros(size(comb_tasks,1),proj_params.dim_manifold);
    vaf_CCA{ds}.vaf_CC2 = zeros(size(comb_tasks,1),proj_params.dim_manifold);
    vaf_CCA{ds}.norm_vaf_CC1 = zeros(size(comb_tasks,1),proj_params.dim_manifold);
    vaf_CCA{ds}.norm_vaf_CC2 = zeros(size(comb_tasks,1),proj_params.dim_manifold);
    
    
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

        % I'm still trying to figure this one out, but I got something from
        % the Statistica website: From:
        % http://www.statsoft.com/Textbook/Canonical-Analysis, which seems
        % to match another two sebsites:
        % http://core.ecu.edu/psyc/wuenschk/MV/Canonical/Canonical.docx and
        % http://sunsite.univie.ac.at/textbooks/statistics/stcanan.html 

        
        % It seems that the average of the square loadings gives you an
        % estimate of VAF, although I'm not sure I got their funky
        % nomenclature right

        vaf1 = zeros(1,proj_params.dim_manifold);
        vaf2 = zeros(1,proj_params.dim_manifold);

        for k = 1:proj_params.dim_manifold

            % This is what they propose ...
            vaf1(k) = mean(A(:,k).^2);
            vaf2(k) = mean(B(:,k).^2);

        %     % But what about this?
        %     vaf1(k) = sum(A(:,k).^2);
        %     vaf2(k) = sum(B(:,k).^2);
        end

        % compute percentage of variance explained
        vaf1 = vaf1/sum(vaf1);
        vaf2 = vaf2/sum(vaf2);

        
        % but this is relative to the scores, but not all the scores reflect the
        % same VAF in the neural data... multiply them both?
        vaf_manifold1 = svd(cov(ssc1))/sum(svd(cov(ssc1)));
        vaf_manifold2 = svd(cov(ssc2))/sum(svd(cov(ssc2)));

        vaf_CCs1 = (vaf1.*vaf_manifold1')/sum(vaf1.*vaf_manifold1');
        vaf_CCs2 = (vaf2.*vaf_manifold2')/sum(vaf2.*vaf_manifold2');
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store the results
        
        vaf_CCA{ds}.vaf_CC1(c,:) = vaf1;
        vaf_CCA{ds}.vaf_CC2(c,:) = vaf2;
        vaf_CCA{ds}.norm_vaf_CC1(c,:) = vaf_CCs1;
        vaf_CCA{ds}.norm_vaf_CC2(c,:) = vaf_CCs2;    
    
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

all_norm_vaf_CC1 = cell2mat( cellfun(@(x) x.norm_vaf_CC1, vaf_CCA, 'UniformOutput',false)' );
all_norm_vaf_CC2 = cell2mat( cellfun(@(x) x.norm_vaf_CC2, vaf_CCA, 'UniformOutput',false)' );

% pool all comparisons together
all_vaf_CC = [all_vaf_CC1; all_vaf_CC2];
all_norm_vaf_CC = [all_norm_vaf_CC1; all_norm_vaf_CC2];


% -------------------------------------------------------------------------
% 2. Wrist tasks and Reach-to-grasp tasks separately

% define vars
all_vaf_CC1_wrist = [];
all_vaf_CC2_wrist = [];
all_norm_vaf_CC1_wrist = [];
all_norm_vaf_CC2_wrist = [];
all_vaf_CC_wrist = [];
all_norm_vaf_CC_wrist = [];

all_vaf_CC1_reach = [];
all_vaf_CC2_reach = [];
all_norm_vaf_CC1_reach = [];
all_norm_vaf_CC2_reach = [];
all_vaf_CC_reach = [];
all_norm_vaf_CC_reach = [];


for d = 1:length(wrist_ds)
    
    all_vaf_CC1_wrist = [all_vaf_CC1_wrist; vaf_CCA{wrist_ds(d)}.vaf_CC1]; %#ok<*AGROW>
    all_vaf_CC2_wrist = [all_vaf_CC2_wrist; vaf_CCA{wrist_ds(d)}.vaf_CC2];
    
    all_norm_vaf_CC1_wrist = [all_norm_vaf_CC1_wrist; vaf_CCA{wrist_ds(d)}.norm_vaf_CC1];
    all_norm_vaf_CC2_wrist = [all_norm_vaf_CC2_wrist; vaf_CCA{wrist_ds(d)}.norm_vaf_CC2];
end


for d = 1:length(reach_ds)
    
    all_vaf_CC1_reach = [all_vaf_CC1_reach; vaf_CCA{reach_ds(d)}.vaf_CC1];
    all_vaf_CC2_reach = [all_vaf_CC2_reach; vaf_CCA{reach_ds(d)}.vaf_CC2];
    
    all_norm_vaf_CC1_reach = [all_norm_vaf_CC1_reach; vaf_CCA{reach_ds(d)}.norm_vaf_CC1];
    all_norm_vaf_CC2_reach = [all_norm_vaf_CC2_reach; vaf_CCA{reach_ds(d)}.norm_vaf_CC2];
end


% pool all comparisons together
all_vaf_CC_wrist = [all_vaf_CC1_wrist; all_vaf_CC2_wrist];
all_norm_vaf_CC_wrist = [all_norm_vaf_CC1_wrist; all_norm_vaf_CC2_wrist];

all_vaf_CC_reach = [all_vaf_CC1_reach; all_vaf_CC2_reach];
all_norm_vaf_CC_reach = [all_norm_vaf_CC1_reach; all_norm_vaf_CC2_reach];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT


% -------------------------------------------------------------------------
% BOTH SETS OF TASKS (WRIST AND REACH-TO-GRASP) TOGETHER

% VAF of the CCs
figure, hold on
for c = 1:size(all_vaf_CC1,1)
    plot(all_vaf_CC1(c,:)*100,'color',[.7 .7 .7])
    plot(all_vaf_CC2(c,:)*100,'color',[.7 .7 .7])
end
plot(mean(all_vaf_CC,1)*100,'k','linewidth',2)
plot(mean(all_vaf_CC,1)*100+std(all_vaf_CC,0,1)*100,'-.k','linewidth',2)
plot(mean(all_vaf_CC,1)*100-std(all_vaf_CC,0,1)*100,'-.k','linewidth',2)
set(gca,'TickDir','out','FontSize',14), box off
ylabel('Variance accounted for')
xlabel('Neural mode')


% VAF in the manifold
figure, hold on
for c = 1:size(all_vaf_CC1,1)
    plot(all_norm_vaf_CC1(c,:)*100,'color',[.7 .7 .7])
    plot(all_norm_vaf_CC2(c,:)*100,'color',[.7 .7 .7])
end
% plot(mean(all_norm_vaf_CC,1),'k','linewidth',2)
% plot(mean(all_norm_vaf_CC,1)+std(all_norm_vaf_CC,0,1),'-.k','linewidth',2)
% plot(mean(all_norm_vaf_CC,1)-std(all_norm_vaf_CC,0,1),'-.k','linewidth',2)
errorbar(1:proj_params.dim_manifold,mean(all_norm_vaf_CC,1)*100,std(all_norm_vaf_CC,0,1)*100,...
    'linestyle','none','linewidth',2,'color','k','marker','.','markersize',32)
set(gca,'TickDir','out','FontSize',14), box off
ylabel('Normalized variance accounted for (%)')
xlabel('Neural mode')
    


% -------------------------------------------------------------------------
% WRIST TASKS AND REACH-TO-GRASP TASKS SEPARATELY


% VAF of the CCs
% figure, hold on
% for c = 1:size(all_vaf_CC1,1)
%     plot(all_vaf_CC1_wrist(c,:),'color',[.7 .7 .7])
%     plot(all_vaf_CC2_wrist(c,:),'color',[.7 .7 .7])
% end
% errorbar(1:proj_params.dim_manifold,mean(all_vaf_CC,1),std(all_norm_vaf_CC,0,1),...
%     'linestyle','none','linewidth',2,'color','k','marker','.','markersize',32)
% set(gca,'TickDir','out','FontSize',14), box off
% ylabel('Variance accounted for')
% xlabel('Neural mode')


% VAF in the manifold
figure, hold on
plot(all_norm_vaf_CC1_wrist(1,:)*100,'color',[.7 .7 .7])
plot(all_norm_vaf_CC1_reach(1,:)*100,'color',[1 .7 .3])
for c = 2:size(all_vaf_CC1_wrist,1)
    plot(all_norm_vaf_CC2_wrist(c,:)*100,'color',[.7 .7 .7])
end
for c = 2:size(all_vaf_CC_reach,1)
    plot(all_norm_vaf_CC_reach(c,:)*100,'color',[1 .7 .3])
end
errorbar(1:proj_params.dim_manifold,mean(all_norm_vaf_CC_wrist,1)*100,std(all_norm_vaf_CC_wrist,0,1)*100,...
    'linestyle','none','linewidth',2,'color','k','marker','.','markersize',32)
errorbar(1:proj_params.dim_manifold,mean(all_norm_vaf_CC_reach,1)*100,std(all_norm_vaf_CC_reach,0,1)*100,...
    'linestyle','none','linewidth',2,'color','r','marker','.','markersize',32)
set(gca,'TickDir','out','FontSize',14), box off
ylabel('Normalized variance accounted for (%)')
xlabel('Neural mode')
legend('wrist','reach-to-grasp','Location','NorthEast'), legend boxoff
    
