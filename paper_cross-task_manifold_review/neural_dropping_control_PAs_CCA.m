%
% Neural dropping controls, done properly
%
% 1. Compute PAs after dropping a certain % of channels (i.e. making them
% zero). Compare the results to the confidence bounds obtained with our
% shuffling control (remember it was virtually identical to TME, and this
% is just more convenient) with all units
%
% 2. Compute CCs using the same approach. Significance threshold computing
% by shuffling over time + smoothing, as in the paper
%

clearvars -except datasets
close all;


% Datasets and task per dataset --these are the 4 examples in Suppl Fig
% 3e,f
ds = [1 8 11 4];
tasks = [3 4 2 2];

% manifold dimensionality
mani_dims = 12;

% percentages of channels to drop
perc_drop = [.5 .4 .3 .2 .1];
% repetitions for each % drop
reps = 100; 

% number of shuffles bootstrapping
n_shuff_CCA = 5000; % as in the paper


% load all the data
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% Load significance PA thresholds obtained after re-shuffling. Remember this
% was almost equal to the TME one
% -- Also, these thresholds were computed for all units
if ~exist('angle_non_orth','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets_original_submission.mat');
end

% bin_size and kernel SD for smoothing the shuffled LVs (used sd = 50 ms, in
% the paper)
bin_size = datasets{ds(1)}.stdata{1}.target{1}.bin_size;
kernel_sd = 0.05;
P_th_CCs = 0.001;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PAs after dropping units
%
% The idea here is to turn off random subsets of channels, so as to keep
% the space dimensionality constant
%

PAs_shuf = zeros(reps,length(perc_drop),mani_dims,length(ds));

for d = 1:length(ds)
   
    fr = datasets{ds(d)}.stdata{tasks(d)}.target{end}.neural_data.conc_smoothed_fr;
    
    n_units = size(fr,2);
    
    for p = 1:length(perc_drop)
        
        for r = 1:reps

            % drop a certain percentage
            rand1 = randperm(n_units);
            drop1 = rand1(1:floor(length(rand1)*perc_drop(p)));

            rand2 = randperm(n_units);
            drop2 = rand2(1:floor(length(rand2)*perc_drop(p)));

            fr1 = fr;
            fr2 = fr;

            fr1(:,drop1) = zeros(size(fr1,1),length(drop1));
            fr2(:,drop2) = zeros(size(fr2,1),length(drop2));

            w1 = pca(fr1);
            w2 = pca(fr2);

            w1 = w1(:,1:mani_dims);
            w2 = w2(:,1:mani_dims);

            PAs_shuf(r,p,:,d) = principal_angles(w1,w2);
        end
    end
end


% stats for the plot
mn_PA = mean(PAs_shuf,1);
sd_PA = std(PAs_shuf,0,1);


cols = parula(length(perc_drop));


figure,
for d = 1:length(ds)
    
    ano = angle_non_orth(:,1,ds(d));
    
    subplot(1,length(ds),d), hold on
    for p = 1:length(perc_drop) 
        plot(rad2deg(squeeze(mn_PA(:,p,:,d))),'color',cols(p,:),'linewidth',1.5)
    end
    plot(ano,'color',[.65 .65 .65],'linewidth',1.5,'linestyle','-.')
    set(gcf,'color','w')
    set(gca,'FontSize',14,'TickDir','out')
    title([datasets{ds(d)}.monkey ' - ' datasets{ds(d)}.labels{tasks(d)} ' - N = ' num2str(length(datasets{ds(d)}.neural_chs))]);
    
    if d == 1
        for p = 1:length(perc_drop)
            lgnd{p} = [num2str(perc_drop(p)*100) ' % drop'];
        end
        lgnd{length(lgnd)+1} = 'P<0.001';
        
        ylabel('Principal angle (deg)');
    end
    legend(lgnd,'Location','NorthWest'),legend boxoff
    xlabel('Neural mode'), xlim([0 mani_dims]), ylim([0 90])
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CCs after dropping units
%
% The idea here is to turn off random subsets of channels, so as to keep
% the space dimensionality constant
%

CCs_shuf = zeros(reps,length(perc_drop),mani_dims,length(ds));
CC_th = zeros(length(ds),mani_dims);

for d = 1:length(ds)
   
    fr = datasets{ds(d)}.stdata{tasks(d)}.target{end}.neural_data.conc_smoothed_fr;
    
    n_units = size(fr,2);
    
    % ---------------------------------------------------------------------
    % Dropout analysis
    
    for p = 1:length(perc_drop)
        
        for r = 1:reps

            % drop a certain percentage
            rand1 = randperm(n_units);
            drop1 = rand1(1:floor(length(rand1)*perc_drop(p)));

            rand2 = randperm(n_units);
            drop2 = rand2(1:floor(length(rand2)*perc_drop(p)));

            fr1 = fr;
            fr2 = fr;

            fr1(:,drop1) = zeros(size(fr1,1),length(drop1));
            fr2(:,drop2) = zeros(size(fr2,1),length(drop2));

            [~, sc1] = pca(fr1);
            [~, sc2] = pca(fr2);

            sc1 = sc1(:,1:mani_dims);
            sc2 = sc2(:,1:mani_dims);

            [~,~, CCs_shuf(r,p,:,d)] = canoncorr(sc1,sc2);
        end
    end
    
    
    % ---------------------------------------------------------------------
    % Do shuffling control
    
    all_CCs_boots = zeros(n_shuff_CCA,mani_dims);
     
    for s = 1:n_shuff_CCA
       
        sc = datasets{ds(d)}.stdata{tasks(d)}.target{end}.neural_data.dim_red.scores;
        sc = sc(:,1:mani_dims);
        
        rnd_sc = sc(randperm(size(sc,1)),:);

        % smooth these shuffled activity patterns to have
        % smoothed random trajectories
        smooth_rnd_sc = smooth_data( rnd_sc , bin_size, kernel_sd );

        [~,~,all_CCs_boots(s,:)] = canoncorr(smooth_rnd_sc,sc);
    end
    
    CC_th(d,:) = prctile(all_CCs_boots,(1-P_th_CCs)*100);
end


% stats for the plot
mn_CC = mean(CCs_shuf,1);
sd_CC = std(CCs_shuf,0,1);



figure,
for d = 1:length(ds)
    
    ano = angle_non_orth(:,1,ds(d));
    
    subplot(1,length(ds),d), hold on
    for p = 1:length(perc_drop) 
        plot(squeeze(mn_CC(:,p,:,d)),'color',cols(p,:),'linewidth',1.5)
    end
    plot(CC_th(d,:),'color',[.65 .65 .65],'linewidth',1.5,'linestyle','-.')
    set(gcf,'color','w')
    set(gca,'FontSize',14,'TickDir','out')
    title([datasets{ds(d)}.monkey ' - ' datasets{ds(d)}.labels{tasks(d)} ' - N = ' num2str(length(datasets{ds(d)}.neural_chs))]);
    
    if d == 1
        for p = 1:length(perc_drop)
            lgnd{p} = [num2str(perc_drop(p)*100) ' % drop'];
        end
        lgnd{length(lgnd)+1} = 'P<0.001';
        ylabel('CC latent activity')
    end
    legend(lgnd,'Location','West'),legend boxoff
    xlabel('Neural mode'), xlim([0 mani_dims]), ylim([0 1])
end