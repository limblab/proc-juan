%
% New control for PA where we shuffle the weights of the linear
% combinations
% 


% gcp;

% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% Load params for CCA --- these are the same ones used in the original
% submission of the paper
proj_params = batch_compare_manifold_projs_defaults();

% Plot for each comparison? 
plot_per_comp_flg = false;

% -------------------------------------------------------------------------
% Parameters for the control analysis

% Choose type of control: 
%  - 'shuffle_weights';
%  - 'shuffle_across_dimensions_and_targets'
%  - 'shuffle_across_neurons_and_targets' -> same as the previous one,
%  but shuffles the neural activity rather than the mode dynamics
%  - 'shuffleSVD_U'
%  - 'shuffle_over_time' (not fully implemented)  
control = 'shuffle_across_neurons_and_targets'; 

% Number of shuffles
n_shuffles = 1000;
P_th = 0.001;


% Load our shuffled principal angles
load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets_original_submission.mat')


% -------------------------------------------------------------------------
% Parameters for analysis

% % Make all the analysis windows 700 ms
% win_length = 0.7;
% for i = 1:size(proj_params.time_win,1)
%     proj_params.time_win(i,2) = proj_params.time_win(i,1)+win_length;
% end

% overwrite manifold dimension, if necessary
proj_params.dim_manifold = 12;


% -------------------------------------------------------------------------

% define which datasets are wrist and which ones are reach-to-grasp,
% separately
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define matrices to store the results
signif_PA.all_actual_PAs = [];
signif_PA.wrist_flg = [];
signif_PA.session_nbr = [];
signif_PA.signif_th = [];
signif_PA.our_shuffled_PA = [];
if strcmp(control,'shuffle_across_dimensions_and_targets')
    signif_PA.signif_th1 = [];
    signif_PA.signif_th2 = [];
end
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

    
    % Do for each task comparison
    for c = 1:size(comb_tasks,1)
    
        
        disp(['Generating random distributions for dataset ' num2str(ds) ...
            ' task comparison ' num2str(c) '/' num2str(size(comb_tasks,1))]);
        
        
        % idx tasks
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
                
        
        n_units = length(stdata{t1}.target{1}.neural_data.neural_chs);

        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE PAs OF THE ACTUAL DATA

        actualPAs = principal_angles(datasets{ds}.dim_red_FR{t1}.w(:,1:proj_params.dim_manifold),...
                                        datasets{ds}.dim_red_FR{t2}.w(:,1:proj_params.dim_manifold));



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
        all_PA_shuff = zeros(n_shuffles,proj_params.dim_manifold);
        
        
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
                    rnd_w1 = w1(idx_shuffle);
                    rnd_w2 = w2(idx_shuffle);
                    
                    rnd_w1 = rnd_w1(:,1:proj_params.dim_manifold);
                    rnd_w2 = rnd_w2(:,1:proj_params.dim_manifold);
                    
                    % Compute the CC
                    all_PA_shuff(s,:) = principal_angles(rnd_w1,rnd_w2);
                end
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONTROL 2 - SHUFFLE OVER TIME, AND THEN SMOOTH TO HAVE SIMILAR
        % DYNAMICS
        
            case 'shuffle_across_dimensions_and_targets'
                
                for s = 1:n_shuffles
                                      
%                     % shuffle the scores across dimensions and trials
%                     % psc1 has size time x trials x neurons matrix
%                     psc1_manifold = psc1(:,:,1:proj_params.dim_manifold);
%                     
%                     rpsc1_manifold = reshape(psc1_manifold,size(psc1_manifold,1),[]);                    
%                     shuffled_rpsc1_manifold = rpsc1_manifold(:,randperm(size(psc1_manifold,2)*size(psc1_manifold,3)));
%                     rnd_ssc1 = reshape(shuffled_rpsc1_manifold,size(psc1_manifold,1),size(psc1_manifold,2),size(psc1_manifold,3));
%                     rnd_ssc1 = reshape(rnd_ssc1,size(psc1_manifold,1)*size(psc1_manifold,2),[]);
%                       
%                     % do the same for the other task
%                     psc2_manifold = psc2(:,:,1:proj_params.dim_manifold);
%                     
%                     rpsc2_manifold = reshape(psc2_manifold,size(psc2_manifold,1),[]);                    
%                     shuffled_rpsc2_manifold = rpsc2_manifold(:,randperm(size(psc2_manifold,2)*size(psc2_manifold,3)));
%                     rnd_ssc2 = reshape(shuffled_rpsc2_manifold,size(psc2_manifold,1),size(psc2_manifold,2),size(psc2_manifold,3));
%                     rnd_ssc2 = reshape(rnd_ssc2,size(psc2_manifold,1)*size(psc2_manifold,2),[]);
%                     
%                     
%                     [~,~,all_PA_shuff(s,:)] = canoncorr(rnd_ssc1,rnd_ssc2);
                end
               
                
            case 'shuffle_across_neurons_and_targets'

                for s = 1:n_shuffles

                    % shuffle the scores across dimensions and trials
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
                    rnd_w1 = pca(rnd_rpfr1);
                    rnd_w2 = pca(rnd_rpfr2);

                    rnd_w1 = rnd_w1(:,1:proj_params.dim_manifold);
                    rnd_w2 = rnd_w2(:,1:proj_params.dim_manifold);

                    all_PA_shuff(s,:) = principal_angles(rnd_w1,rnd_w2);
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
%                     all_PA_shuff(s,:) = calc_r(rnd_ssc1,rnd_ssc2);

                end                

        end
        
        
        % -----------------------------------------------------------------
        % Store results
        
        signif_PA.all_actual_PAs = [signif_PA.all_actual_PAs; actualPAs];
        if ismember(ds,wrist_ds), wyn = 1; else wyn = 0; end
        signif_PA.wrist_flg = [signif_PA.wrist_flg; wyn];
        signif_PA.session_nbr = [signif_PA.session_nbr; ds];
        
        
        % Our shuffled control
        our_shuffled_PA = angle_non_orth(:,1,space_dim==n_units)';
        signif_PA.our_shuffled_PA = our_shuffled_PA;
        
        t_signif_th = prctile(all_PA_shuff,P_th*100);
        signif_PA.signif_th = [signif_PA.signif_th; t_signif_th]; %#ok<*AGROW>
        
        ctr = ctr + 1;
        
        % -----------------------------------------------------------------
        % plot
        
        if plot_per_comp_flg
            figure,hold on %#ok<UNRCH>
            plot(rad2deg(actualPAs),'k','linewidth',2)
            p1 = plot(rad2deg(t_signif_th),'color',[.5 0 .5],'linewidth',2.5);
            p2 = plot(our_shuffled_PA,'color',[.6 .6 .6],'linewidth',2.5,'linestyle','-.');
            plot(rad2deg(all_PA_shuff)','color',[216 190 216]/255)
            ylim([0 90]),xlim([0 proj_params.dim_manifold]);
            set(gca,'TickDir','out','FontSize',14), box off
            xlabel('Neural mode'),ylabel('Principal angle (deg)')
            legend('actual','rnd PSTHs','our control','surrogates rnd PSTHs','Location','NorthWest'),legend boxoff
            uistack(p1,'top');uistack(p2,'top');
            
            pause; close
        end
    end
end



% -------------------------------------------------------------------------
%% COMPARISON OF OUR SHUFFLING CONTROL AND THIS CONTROL



lfit = polyfit(rad2deg(signif_PA.signif_th),our_shuffle_th,1);
xfit4plot = [rad2deg(min(min(min(signif_PA.signif_th),min(our_shuffle_th))))-5 90];
yfit4plot = polyval(lfit,xfit4plot);

% compute correlation
[r, Pr] = corr(reshape(rad2deg(signif_PA.signif_th),[],1),reshape(our_shuffle_th,[],1));

hf = figure; hold on
plot([0 90],[0 90],'color',[.6 .6 .6])
plot(xfit4plot,yfit4plot,'k','linewidth',1.5)
plot(rad2deg(signif_PA.signif_th),our_shuffle_th,'.k','markersize',12)
legend('perfectly equal','method match','Location','SouthEast'),legend boxoff
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(hf, 'color', [1 1 1]);
text(10,75,[num2str(lfit(2),3) '+' num2str(lfit(1),3) '·x'],'FontSize',14)
text(10,65,['r=' num2str(r,3) '; P=' num2str(Pr,3)],'FontSize',14)
xlabel('New control'); ylabel('Our random sampling threshold')
xlim([0 90]),ylim([0 90])

