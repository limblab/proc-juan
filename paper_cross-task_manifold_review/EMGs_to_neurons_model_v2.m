%
% Simple model that generates neural activity assuming it is a weighted
% linear combination of the EMGs, and then compares the simulated activity
% between two tasks (using a fixed model) with CCA
%
% Much of the model stuff inspired by Perich & Miller EBR 2017
%


clc;

mani_dim = 12;
ds_to_use = [1:3 7:9];


add_noise = true; % add white noise to the Poisson process with neuron-dependent variable 'gain'
add_measure_noise = false; % add Gaussian noise to the simulated spikes to simulate measurement noise
do_smoothing = true; % make spikes really smooth --not necessary because spikes are simulated from smoothed EMGs
scale_frs = true; % scale the FRs so the mean for each simulated neuron matches the real data
add_lag = true;

lag_range = [50 150];

% gain for the WN applied to the lambdas
noise_amplif = .1;

% params for Gaussian smoothing
bin_size = 0.02;
kernel_sd = 0.05;

% one plot per comparison?
comp_plots = false;


% Define vars to store main results
all_model_CCs = [];
all_real_CCs = [];
all_model_screes = [];
all_real_screes = [];


% -------------------------------------------------------------------------

% The real deal: do!
for d = 1:length(ds_to_use)
    
    % get dataset idx
    t_ds = ds_to_use(d);
    
    comb_tasks = nchoosek(1:length(datasets{t_ds}.labels),2);
    
    % -------------------------------------------------------------------------
    % Prepare data for CCA

    % a) equalize trial duration across all tasks
    stdata   = equalize_single_trial_dur( datasets{t_ds}.stdata );
    % b) equalize number of trials for all targets of a given task
    for i = 1:length(stdata)
        stdata{i} = equalize_nbr_trials_p_target( stdata{i} );
    end
    % c) equalize number of trials across tasks
    stdata = equalize_nbr_trials_across_tasks( stdata );


    % -------------------------------------------------------------------------
    % Model Parameters (most adopted from Matt's paper)
    %
    % --the same for all the tasks we're comparing

    % For each neuron, the model is: lambda = w_0 + Sum_i{w_i*alpha_i}
    
    % same number of neurons and EMGs as in the data
    n_neurons = length(datasets{t_ds}.neural_chs);
    n_emgs = length(datasets{t_ds}.chosen_emgs);
    
    % ---------------------------------------------------------------------
    % Generate weights
    
    % w_0: drawn from a Gaussian distribution. Have to range from 0 to 0.1
    w_0 = abs(randn(n_neurons,1));
    
    % Scale is so it ranges from 0 to 0.1
    w_0 = w_0/20; % 95% of the weights should now range be within [0, 0.1]
    w_0(w_0>.1) = .1; % make weights > 0.1 be 0.1
    
    % w_1: drawn from a Gaussian. Have to go from -1 to 1
    w_1 = randn(n_neurons,n_emgs);
    
    % scale them so they range from -1 to 1
    w_1 = w_1/2; % 95% of the weights should now range be within [-1, 1]
    % Trick: scale the weights in w_i so they go from -1 to 1
    w_1(w_1>1) = 1;
    w_1(w_1<-1) = -1;

    
    % Upper and lower bounds for lambda 
    up_th_lambda = 1; % upper limit for lambda -- we don't want it too large
    low_th_lambda = 0; % lower limit for lambda -- we don't want it to be negative!!

    
    % Generate noise gains, if we add white noise to the lambdas
    noise_gain = noise_amplif*randn(n_neurons,1)/2;
    
    % Generate gains for measurement noise
    measure_noise_gain = randn(n_neurons,1)/20;
    
    % Generate lags
%     % Normally distributed
%     r_emg_neuron_lag = randn(n_neurons,1);
%     emg_neuron_lag = lag_range(1) + diff(lag_range)/2 + ( diff(lag_range)/2 * r_emg_neuron_lag/3 );
    % Uniformly distributed    
    r_emg_neuron_lag = rand(n_neurons,1);
    emg_neuron_lag = lag_range(1) + r_emg_neuron_lag * diff(lag_range);
    emg_neuron_lag = emg_neuron_lag / 1000;
    lag_bins = round( emg_neuron_lag / bin_size );
    
% -------------------------------------------------------------------------
% SIMULATION, one per task comparison

    for c = 1:size(comb_tasks,1)
    
        % Task idx
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
        
        % Get the data
        real_neurons1 = stdata{t1}.target{end}.neural_data.conc_smoothed_fr';
        real_emg1 = stdata{t1}.target{end}.emg_data.conc_emg';
        real_neurons2 = stdata{t2}.target{end}.neural_data.conc_smoothed_fr';
        real_emg2 = stdata{t2}.target{end}.emg_data.conc_emg';

        n_bins = size(real_neurons1,2);
        
        % 1a. Compute the lambdas -without lag
        if ~add_lag
            lambda1 = w_0 + w_1*real_emg1;
            lambda2 = w_0 + w_1*real_emg2;
        end
        % 1b. Compute the lambdas with lag
        if add_lag
            
            lambda1 = zeros(n_neurons,n_bins);
            lambda2 = zeros(n_neurons,n_bins);
            
            for n = 1:n_neurons
                idx_beg = lag_bins(n)+1;
                idx_end = n_bins-lag_bins(n);
                
                lambda1(n,1:idx_end) = w_0(n) + w_1(n,:)*real_emg1(:,idx_beg:end);
                lambda2(n,1:idx_end) = w_0(n) + w_1(n,:)*real_emg2(:,idx_beg:end);
            end
            
            % trim the last spikes that will be zeroes because of the lag
            lambda1(:,end-round(lag_range(2)/bin_size/1000)+1:end) = [];
            lambda2(:,end-round(lag_range(2)/bin_size/1000)+1:end) = [];
        end

        % 2 Add noise?
        if add_noise
           lambda1 = lambda1 + noise_gain.*randn(size(lambda1)); 
           lambda2 = lambda2 + noise_gain.*randn(size(lambda2)); 
        end
        
        % 3. Apply upper and lower limits
        lambda1 = abs(lambda1);
        lambda2 = abs(lambda2);
        
        lambda1(lambda1>up_th_lambda) = up_th_lambda;
        lambda1(lambda1<low_th_lambda) = low_th_lambda;
        lambda2(lambda2>up_th_lambda) = up_th_lambda;
        lambda2(lambda2<low_th_lambda) = low_th_lambda;
        
        
        % 4. Compute the firing rates
        fr1 = poissrnd(lambda1);
        fr2 = poissrnd(lambda2);

       
        % 5a. Preprocess them
        % square root transform
        fr1 = sqrt(fr1);
        fr2 = sqrt(fr2);
        
%         % 5b. Add measurement noise? -- THIS USED TO BE BEFORE THE SQRT WHICH WAS WRONG
%         if add_measure_noise
%             fr1 = fr1 + measure_noise_gain.*randn(size(fr1)); % THIS EQ IS WRONG
%             fr2 = fr2 + measure_noise_gain.*randn(size(fr2));
%         end
        
        
        % 6. Gaussian smoothing
        if do_smoothing
            fr1 = smooth_data(fr1',bin_size,kernel_sd)';
            fr2 = smooth_data(fr2',bin_size,kernel_sd)';
        end
        
        
        % 7. Scale the firing rates, so their mean matches the real data
        if scale_frs
            
            gain1 = mean(real_neurons1,2)./mean(fr1,2);
            gain2 = mean(real_neurons2,2)./mean(fr2,2);
            % trick -- if gain is Inf make it 1 ;-)
            gain1(isinf(gain1)) = 1;
            gain2(isinf(gain2)) = 1;
            % apply the scaling
            scaled_fr1 = gain1.*fr1;
            scaled_fr2 = gain2.*fr2;
        else
            scaled_fr1 = fr1;
            scaled_fr2 = fr2;
        end


        % -----------------------------------------------------------------
        % Do PCA of the simulated FRs
        [lv1,W1,eigenvals1] = pca(scaled_fr1);
        [lv2,W2,eigenvals2] = pca(scaled_fr2);

        % Normalized cumulative vaf
        scree1 = cumsum(eigenvals1)/sum(eigenvals1);
        scree2 = cumsum(eigenvals2)/sum(eigenvals2);
        % keep the max between 20 and mani_dim
        scree1 = scree1(1:max(20,mani_dim))';
        scree2 = scree2(1:max(20,mani_dim))';

        % -----------------------------------------------------------------
        % CCA to compare the dynamics
        
        % Keep dimensions in latent activity that we want
        lv1 = lv1(:,1:mani_dim);
        lv2 = lv2(:,1:mani_dim);
        
        % CCA between the simulated neural data
        [~,~,CCs] = canoncorr(lv1,lv2);

        % CCA between the real data --note that this is not the same windows as for
        % the real paper
        [~,~,CCs_real] = canoncorr(stdata{t1}.target{end}.neural_data.dim_red.scores(:,1:mani_dim), ...
                                    stdata{t2}.target{end}.neural_data.dim_red.scores(:,1:mani_dim));


        % -----------------------------------------------------------------
        % Store results
        
        all_model_CCs = [all_model_CCs; CCs]; %#ok<*AGROW>
        all_real_CCs = [all_real_CCs; CCs_real];
        all_model_screes = [all_model_screes; scree1; scree2];
            
        % -----------------------------------------------------------------
        % Plots

        if comp_plots 

            % Scree plot PCA simulated data
            figure,hold on,
            plot(cumsum(eigenvals1(1:20))/sum(eigenvals1)*100,'linewidth',1.5,'color','k')
            plot(cumsum(eigenvals2(1:20))/sum(eigenvals2)*100,'linewidth',1.5,'color','r')
            set(gca,'TickDir','out','FontSize',14), box off,xlim([0 20]),ylim([0 100])
            ylabel('Neural variance explained model (%)'),xlabel('Neural modes')
            legend('Task 1','Task 2','Location','SouthEast'), legend boxoff

            % CC for the simulated and neural data
            figure,hold on
            plot(CCs,'linewidth',1.5,'color','k')
            plot(CCs_real,'linewidth',1.5,'color',[.6 .6 .6])
            set(gca,'TickDir','out','FontSize',14), box off, xlim([0 mani_dim]),ylim([0 1])
            ylabel('Canonical correlation dynamics'), xlabel('Neural mode')
            legend('Model','Real','Location','NorthEast'), legend boxoff
        end        
    end
    
    % ---------------------------------------------------------------------
    % Compute the scree plots for the real data, for comparison
        
    for t= 1:length(datasets{t_ds}.labels)
        
        % compute real scree
        real_scree = cumsum(datasets{t_ds}.dim_red_FR{t}.eigen/sum(datasets{t_ds}.dim_red_FR{t}.eigen))';
        real_scree = real_scree(1:max(mani_dim,20));
        
        all_real_screes = [all_real_screes; real_scree];
    end
end


% -------------------------------------------------------------------------
% Summary plot       


% some summary calculations
m_model_CCs = mean(all_model_CCs,1);
sd_model_CCs = std(all_model_CCs,0,1);
m_real_CCs = mean(all_real_CCs,1);
sd_real_CCs = std(all_real_CCs,0,1);

m_scree_model = 100*mean(all_model_screes,1);
sd_scree_model = 100*std(all_model_screes,0,1);
m_scree_real = 100*mean(all_real_screes,1);
sd_scree_real = 100*std(all_real_screes,0,1);


% All scree plots --note that there are two times more than CCs because we
% keep each simulated task's scree plot and we simulate the two to compare
% figure,subplot(121),hold on
% plot(100*all_model_screes(1,:)','color','k','linewidth',1.5);
% plot(100*all_real_screes(1,:)','color',[.6 .6 .6],'linewidth',1.5);
% plot(100*all_model_screes(2:end,:)','color','k','linewidth',1.5);
% plot(100*all_real_screes(2:end,:)','color',[.6 .6 .6],'linewidth',1.5);
% set(gca,'TickDir','out','FontSize',14), box off,xlim([0 20]),ylim([0 100])
% ylabel('Cumulative neural VAF (%)'),xlabel('Neural mode')
% legend('Model','Real','Location','SouthEast'), legend boxoff
            
% subplot(122),hold on
% plot(all_model_CCs(1,:)','linewidth',1.5,'color','k')
% plot(all_real_CCs(1,:)','linewidth',1.5,'color',[.6 .6 .6])
% plot(all_model_CCs(2:end,:)','linewidth',1.5,'color','k')
% plot(all_real_CCs(2:end,:)','linewidth',1.5,'color',[.6 .6 .6])
% set(gca,'TickDir','out','FontSize',14), box off, xlim([0 mani_dim]),ylim([0 1])
% ylabel('CC latent activity'), xlabel('Neural mode')
% legend('Model','Real','Location','NorthEast'), legend boxoff


figure,subplot(121),hold on
errorbar(m_scree_model,sd_scree_model,'.-k','markersize',24,'linewidth',1.5)
errorbar(m_scree_real,sd_scree_real,'.','color',[.6 .6 .6],'markersize',24,'linewidth',1.5,'linestyle','-')
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 20]),ylim([0 100])
ylabel('Cumulative neural VAF (%)'),xlabel('Neural mode')
legend('Model','Real','Location','SouthEast'), legend boxoff

subplot(122),hold on
errorbar(m_model_CCs,sd_model_CCs,'.-k','markersize',24,'linewidth',1.5)
errorbar(m_real_CCs,sd_real_CCs,'.','color',[.6 .6 .6],'markersize',24,'linewidth',1.5,'linestyle','-')
set(gca,'TickDir','out','FontSize',14), box off, xlim([0 mani_dim]),ylim([0 1])
ylabel('CC latent activity'), xlabel('Neural mode')
legend('Model','Real','Location','NorthEast'), legend boxoff

set(gcf, 'color', [1 1 1]);


% Statistical test to compare both distribution of CCs

[~, p] = ttest2(reshape(all_real_CCs,1,[]),reshape(all_model_CCs,1,[]));
disp(['Paired t-test between CC distributions for model and real data: ' num2str(p)]);


% Ratio between model and real data CCs

CCratio = reshape(all_real_CCs,1,[])./reshape(all_model_CCs,1,[]);
mn_CCratio = mean(CCratio);
sd_CCratio = std(CCratio);

disp(['Ratio between real data CCs and model CCs: ' num2str(mn_CCratio) ' +/- ' num2str(sd_CCratio)]);
