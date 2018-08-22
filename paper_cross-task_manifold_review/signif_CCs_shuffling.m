%
% New set of controls for CCA where we: 1) shuffle the neural units'
% weights onto the neural modes, 2) shuffle the neural mode dynamics across
% trials and dimensions, 3) shuffle the neural units' activity patterns
% across units and targets, plus other options not fully implemented
% 
%


gcp;

% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% Load params for CCA --- these are the same ones used in the original
% submission of the paper
proj_params = batch_compare_manifold_projs_defaults();

% Plot for each comparison? 
plot_per_comp_flg = true;

% Manifold dimensionality
mani_dims = 12;

% -------------------------------------------------------------------------
% Parameters for the control analysis

% Choose type of control: 
%  - 'shuffle_weights';
%  - 'shuffle_across_dimensions_and_targets'
%  - 'shuffle_across_neurons_and_targets' -> same as the previous one,
%  but shuffles the neural activity rather than the mode dynamics
%  - 'shuffleSVD_U'
%  - 'shuffle_over_time'
control = 'shuffle_over_time'; % 'shuffle_across_neurons_and_targets'; 
invert_random_trials_flg = true;

% Number of shuffles
n_shuffles = 5000;
P_th = 0.001;


% -------------------------------------------------------------------------
% Parameters for analysis

% % Make all the analysis windows 700 ms
% win_length = 0.7;
% for i = 1:size(proj_params.time_win,1)
%     proj_params.time_win(i,2) = proj_params.time_win(i,1)+win_length;
% end

% overwrite manifold dimension, if necessary
proj_params.dim_manifold = mani_dims;


% -------------------------------------------------------------------------

% define which datasets are wrist and which ones are reach-to-grasp,
% separately
wrist_ds = [1:3 7:9]; % [1:3 7:9];
reach_ds = [4:6 10:11];

ds_to_use = [wrist_ds reach_ds];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define matrices to store the results
signif_CC.all_actual_CCs = [];
signif_CC.wrist_flg = [];
signif_CC.session_nbr = [];
signif_CC.signif_th = [];
if strcmp(control,'shuffle_across_dimensions_and_targets')
    signif_CC.signif_th1 = [];
    signif_CC.signif_th2 = [];
end
ctr = 1;


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


    
    % Do for each task comparison
    for c = 1:size(comb_tasks,1)
    
        
        disp(['Generating random distributions for dataset ' num2str(ds) ...
            ' task comparison ' num2str(c) '/' num2str(size(comb_tasks,1))]);
        
        
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
        fr1 = stdata{comb_tasks(c,1)}.target{end}.neural_data.smoothed_fr;
        fr2 = stdata{comb_tasks(c,2)}.target{end}.neural_data.smoothed_fr;

        % turn fri into (time x trials) x neurons matrix
        pfr1 = permute(fr1,[1 3 2]);
        sfr1 = reshape(pfr1,[],size(pfr1,3));

        pfr2 = permute(fr2,[1 3 2]);
        sfr2 = reshape(pfr2,[],size(pfr2,3));
        
        
        % variable to store the CCA after projecting onto the shuffled
        % manifolds
        all_rshuff = zeros(n_shuffles,proj_params.dim_manifold);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONTROL 1 - SHUFFLE WEIGHTS IN PCA MATRIX
        
        switch control
            
            case 'shuffle_weights'
                
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
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONTROL 2 - SHUFFLE OVER TIME, AND THEN SMOOTH TO HAVE SIMILAR
        % DYNAMICS
        
            case 'shuffle_across_dimensions_and_targets'
                
                for s = 1:n_shuffles
                                      
                    % shuffle the scores across dimensions and trials
                    % psc1 has size time x trials x neurons matrix
                    psc1_manifold = psc1(:,:,1:proj_params.dim_manifold);
                    
                    rpsc1_manifold = reshape(psc1_manifold,size(psc1_manifold,1),[]);                    
                    shuffled_rpsc1_manifold = rpsc1_manifold(:,randperm(size(psc1_manifold,2)*size(psc1_manifold,3)));
                    rnd_ssc1 = reshape(shuffled_rpsc1_manifold,size(psc1_manifold,1),size(psc1_manifold,2),size(psc1_manifold,3));
                    rnd_ssc1 = reshape(rnd_ssc1,size(psc1_manifold,1)*size(psc1_manifold,2),[]);
                      
                    % do the same for the other task
                    psc2_manifold = psc2(:,:,1:proj_params.dim_manifold);
                    
                    rpsc2_manifold = reshape(psc2_manifold,size(psc2_manifold,1),[]);                    
                    shuffled_rpsc2_manifold = rpsc2_manifold(:,randperm(size(psc2_manifold,2)*size(psc2_manifold,3)));
                    rnd_ssc2 = reshape(shuffled_rpsc2_manifold,size(psc2_manifold,1),size(psc2_manifold,2),size(psc2_manifold,3));
                    rnd_ssc2 = reshape(rnd_ssc2,size(psc2_manifold,1)*size(psc2_manifold,2),[]);
                    
                    
                    [~,~,all_rshuff(s,:)] = canoncorr(rnd_ssc1,rnd_ssc2);
                end
               
                
            case 'shuffle_across_neurons_and_targets'

                for s = 1:n_shuffles

                    % "invert" the activity of some units in some trials
                    if invert_random_trials_flg
                        % random matrices of 1s and 0s: 1 => invert. Different for
                        % the data from each task
                        idx_invert1 = randi([0 1],[size(pfr1,2) size(pfr1,3)]);
%                         idx_invert2 = randi([0 1],[size(pfr1,2) size(pfr1,3)]);
                        for d2 = 1:size(pfr1,2)
                            for d3 = 1:size(pfr1,3)
                                if idx_invert1(d2,d3)
                                    pfr1(:,d2,d3) = pfr1(end:-1:1,d2,d3);
                                end
%                               if idx_invert2(d2,d3)
%                                   pfr2(:,d2,d3) = pfr2(end:-1:1,d2,d3);
%                               end
                            end
                        end
                    end
                    
                    % shuffle the firing rates across dimensions and trials
                    % psc1 has size time x trials x neurons matrix
                    rpfr1 = reshape(pfr1,size(pfr1,1),[]);
                    shuffled_rpfr1 = rpfr1(:,randperm(size(rpfr1,2)*size(rpfr1,3)));
                    rnd_rpfr1 = reshape(shuffled_rpfr1,size(shuffled_rpfr1 ,1),size(shuffled_rpfr1 ,2),size(shuffled_rpfr1 ,3));
                    rnd_rpfr1 = reshape(rnd_rpfr1,size(pfr1,1)*size(pfr1,2),[]);

                    % do the same for the other task
                    rpfr2 = reshape(pfr2,size(pfr2,1),[]);
                    shuffled_rpfr2 = rpfr2(:,randperm(size(rpfr2,2)*size(rpfr2,3)));
                    rnd_rpfr2 = reshape(shuffled_rpfr2,size(shuffled_rpfr2 ,1),size(shuffled_rpfr2 ,2),size(shuffled_rpfr2 ,3));
                    rnd_rpfr2 = reshape(rnd_rpfr2,size(pfr2,1)*size(pfr2,2),[]);

                    % Do PCA on the neural activity keeping the leading mode
                    % dynamics
                    [~,rnd_ssc1] = pca(rnd_rpfr1);
                    [~,rnd_ssc2] = pca(rnd_rpfr2);

                    rnd_ssc1 = rnd_ssc1(:,1:proj_params.dim_manifold);
                    rnd_ssc2 = rnd_ssc2(:,1:proj_params.dim_manifold);

                    [~,~,all_rshuff(s,:)] = canoncorr(rnd_ssc1,rnd_ssc2);
                end
                
            case 'shuffleSVD_U'
                
                for s = 1:n_shuffles

%                     [U1, S1, V1] = svd(sfr1);
%                     [U2, S2, V2] = svd(sfr2);
% 
%                     V1sh = V1(randperm(numel(V1)));
%                     V1sh = reshape(V1,size(S1,2),[]);
%                     rnd_ssc1 = U1*S1*V1sh';
% 
%                     V2sh = V2(randperm(numel(V2)));
%                     V2sh = reshape(V2,size(S2,2),[]);
%                     rnd_ssc2 = U2*S2*V2sh';
%                     
%                     % Plot FFTs to check they have similar smoothness
%                     bin_size = 1/stdata{comb_tasks(c,1)}.target{1}.bin_size;
%                     l_ssc1 = length(ssc1);
%                     nfft = 2^nextpow2(l_ssc1);
%                     Yrnd_ssc1 = fft(rnd_ssc1,nfft)/l_ssc1;
%                     Y_ssc1 = fft(ssc1,nfft)/l_ssc1;
%                     f = bin_size/2*linspace(0,1,nfft/2);
%                     
%                     figure, hold on
%                     plot(f,2*abs(Y_ssc1(1:nfft/2)),'k'),plot(f,2*abs(Yrnd_ssc1(1:nfft/2)),'-.c')
%                     legend('real data','shuffled'), legend boxoff
%                     xlabel('Frequency (Hz)'),ylabel('Amplitude')
%                     set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
%                     
%                     % 
%                     all_rshuff(s,:) = calc_r(rnd_ssc1,rnd_ssc2);

                end                
            case 'shuffle_over_time'
                
                
                for s = 1:n_shuffles
                    
                    % This didn't work very well
                    % shuffle the scores for one task over time
                    % -- this ain't a great idea, we don't keep the
                    % smoothness at all, even after shuffling!!!
                    rnd_ssc1 = ssc1(randperm(size(ssc1,1)),:);

                    % smooth these shuffled activity patterns to have
                    % smoothed random trajectories
                    smooth_rnd_ssc1 = smooth_data( rnd_ssc1, datasets{ds}.stdata{1}.target{1}.bin_size, 0.05);


    %                 % Plot FFTs to check they have similar smoothness
    %                 bin_size = 1/stdata{comb_tasks(c,1)}.target{1}.bin_size;
    %                 l_ssc1 = length(ssc1);
    %                 nfft = 2^nextpow2(l_ssc1);
    %                 Yrnd_ssc1 = fft(smooth_rnd_ssc1,nfft)/l_ssc1;
    %                 Y_ssc1 = fft(ssc1,nfft)/l_ssc1;
    %                 f = bin_size/2*linspace(0,1,nfft/2);
    %                 
    %                 figure, hold on
    %                 plot(f,2*abs(Y_ssc1(1:nfft/2)),'k'),plot(f,2*abs(Yrnd_ssc1(1:nfft/2)),'-.c')
    %                 legend('real data','shuffled'), legend boxoff
    %                 xlabel('Frequency (Hz)'),ylabel('Amplitude')
    %                 set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')


                    smooth_rnd_ssc1 = smooth_rnd_ssc1(:,1:proj_params.dim_manifold);

                    [~,~,all_rshuff(s,:)] = canoncorr(smooth_rnd_ssc1,ssc2);
                end
        end
        
        
        % -----------------------------------------------------------------
        % Store results
        
        signif_CC.all_actual_CCs = [signif_CC.all_actual_CCs; r];
        if ismember(ds,wrist_ds), wyn = 1; else wyn = 0; end
        signif_CC.wrist_flg = [signif_CC.wrist_flg; wyn];
        signif_CC.session_nbr = [signif_CC.session_nbr; ds];
        
%         % For the 'shuffle_across_dimensions_and_targets' control we do it
%         % for each task we are comparing
%         if strcmp(control,'shuffle_across_dimensions_and_targets')
%             t_signif_th1 = prctile(all_rshuff1,(1-P_th)*100);
%             t_signif_th2 = prctile(all_rshuff2,(1-P_th)*100);
%             signif_CC.signif_th1 = [signif_CC.signif_th1; t_signif_th1]; %#ok<*AGROW>
%             signif_CC.signif_th2 = [signif_CC.signif_th2; t_signif_th2]; %#ok<*AGROW>
%             % take the significance threshold as the max
%             t_signif_th = max(t_signif_th1,t_signif_th2);
%             
%             signif_CC.all_shuff_CCs1{ctr} = all_rshuff1; %#ok<*SAGROW>
%             signif_CC.all_shuff_CCs2{ctr} = all_rshuff2; %#ok<*SAGROW>
%         else
%             t_signif_th = prctile(all_rshuff1,(1-P_th)*100);
%         end
        t_signif_th = prctile(all_rshuff,(1-P_th)*100);
        signif_CC.signif_th = [signif_CC.signif_th; t_signif_th]; %#ok<*AGROW>
        
        ctr = ctr + 1;
        
        % -----------------------------------------------------------------
        % plot
        
        if plot_per_comp_flg
            figure,hold on %#ok<UNRCH>
            plr = plot(r,'k','linewidth',1.5);
            plot(all_rshuff','color',[216 190 216]/255)
%             if exist('all_rshuff2','var')
%                 plot(all_rshuff','color',[255 215 0]/255)
%             end
            plot(t_signif_th,'color',[.5 0 .5],'linewidth',1.5)
%             if exist('all_rshuff2','var')
%                 plot(t_signif_th2','color',[255 140 0]/255,'linewidth',1.5)
%             end
            uistack(plr,'top')
            ylim([0 1]),xlim([0 proj_params.dim_manifold]);
            set(gca,'TickDir','out','FontSize',14), box off
            xlabel('Neural mode'),ylabel('Canonical correlation')
            legend('actual','surrogates')
            title([datasets{ds_to_use(ds)}.monkey ' - ' datasets{ds_to_use(ds)}.date(1:end-4) ' - ' ...
                datasets{ds_to_use(ds)}.labels{comb_tasks(c,1)} ' vs ' datasets{ds_to_use(ds)}.labels{comb_tasks(c,2)}])
            % pause; close
        end
    end
end


% -------------------------------------------------------------------------
%% Summary calculations and figures


% Find number of modes for each task whose CC is above this significance
% threshold

cc_diff = signif_CC.all_actual_CCs-signif_CC.signif_th;
highest_similar_mode = zeros(1,size(cc_diff,1));

for c = 1:size(cc_diff,1)
    
    t_highest_mode = find(cc_diff(c,:)<0,1)-1;
    if ~isempty(t_highest_mode)
        highest_similar_mode(c) = t_highest_mode;
    else
        highest_similar_mode(c) = proj_params.dim_manifold;
    end
end


% Do histogram of the highest modes
x_hist = 1:proj_params.dim_manifold+1;
y_hist = histcounts(highest_similar_mode,x_hist);


% Plot
figure,
bar(x_hist(1:end-1),y_hist,'FaceColor',[.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Highest mode w. similar dynamics'),ylabel('Counts')



% -------------------------------------------------------------------------
% Number of modes above chance level depending on the task comparison
metadata = batch_get_monkey_task_data(datasets);

highest_similar_mode_per_type = zeros(length(unique(metadata.task_pairs.task_pair_nbr)),...
                                proj_params.dim_manifold);

% task comparisons with zero similar modes
no_similar_modes_ctr = 0;
                            
for task_pair = 1:length(metadata.task_pairs.unique_pairs)
    
    trials_this_pair = find(metadata.task_pairs.task_pair_nbr == task_pair);
    
    for c = 1:length(trials_this_pair)
       
        highest_this = find(cc_diff(trials_this_pair(c),:)<0,1) - 1;
        
        if isempty(highest_this)
            highest_similar_mode_per_type(task_pair,proj_params.dim_manifold) = 1 + ...
                highest_similar_mode_per_type(task_pair,proj_params.dim_manifold);
        elseif highest_this == 0
            % there are no similar modes for this task comparison
            no_similar_modes_ctr = no_similar_modes_ctr + 1;
        else
            highest_similar_mode_per_type(task_pair,highest_this) = 1 + ...
                                highest_similar_mode_per_type(task_pair,highest_this);
        end
    end
end


for p = 1:length(metadata.task_pairs.unique_pairs)
    this_legend{p}         = [metadata.task_pairs.unique_pairs{p}{1} ' vs ' ...
        metadata.task_pairs.unique_pairs{p}{2}];
end



% Summary histogram 
figure,
bar(1:proj_params.dim_manifold,highest_similar_mode_per_type','stacked')
if no_similar_modes_ctr > 1
    hold on
    bar(0,no_similar_modes_ctr,'FaceColor',[.7 .7 .7])
    set(gca,'XTick',0:proj_params.dim_manifold);
    xt{1} = 'No';
    for d = 1:proj_params.dim_manifold
        xt{d+1} = num2str(d);
    end
    set(gca,'XTickLabel',xt)
end
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
xlabel('Highest mode w. similar dynamics'),ylabel('Counts')
legend(this_legend,'Location','NorthWest'), legend boxoff
