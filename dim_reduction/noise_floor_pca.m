%
% Estimate noise floor from neural data using the method presented in
% Machens et al., J Neurosci 2010, extended to continuous neural signals,
% rather than trial separated data
%
% 
% ToDo: NEEDS TO BE FINISHED

function noise_floor = noise_floor_pca( single_trial_data, varargin ) 


% nbr iterations for Machens' method to estimate noise variance in the data
nbr_iter            = 1000;
alpha               = 0.95;

% -------------------------------------------------------------------------
% do PCA on the task-related part of the trial only --this corresponds to
% the concatenated smoothed firing rates of the last target in
% single_trial_data

% first create an array with the smoothed FRs, and add time in the first
% one

aux_sfr             = [single_trial_data.target{end}.conc_t', ...
                        single_trial_data.target{end}.neural_data.conc_smoothed_fr];

dim_red_trial_rel   = dim_reduction( aux_sfr, 'pca' );

nbr_chs             = length(dim_red_trial_rel.eigen);


% -------------------------------------------------------------------------
% do PCA on the trial-averaged firing rates

aux_sfr_trial_rel_avg = [single_trial_data.target{end}.t', ...
                        single_trial_data.target{end}.neural_data.smoothed_fr_mn];

dim_red_trial_rel_avg = dim_reduction( aux_sfr_trial_rel_avg, 'pca' );
                

% -------------------------------------------------------------------------
% Machens' method to estimate the noise floor. Note that it is for
% trial-averaged data

% 1. generate average noise for each neuron by subtracting 2 trials (for all
% targets)

% retrieve number of trials per target
nbr_trials_p_tgt    = cell2mat(cellfun(@(x) size(x.neural_data.fr,3), single_trial_data.target, ...
                        'UniformOutput', false));
nbr_trials_p_tgt    = nbr_trials_p_tgt(1:end-1); % the last one is all concatenated trials

% do n times

for j = 1:nbr_iter
    
    % generate a random pair of numbers which are the trials that will be used
    % to generate the noise estmate
    trial_nbrs    	= datasample(1:min(nbr_trials_p_tgt),2,'Replace',false);

    % fill the firing rate matrices
    sfr_noise    	= [];
    sfr_noise2     	= [];
    for i = 1:length(single_trial_data.target)-1
        sfr_noise   = [sfr_noise; single_trial_data.target{i}.neural_data.smoothed_fr(:,:,trial_nbrs(1))];
        sfr_noise2  = [sfr_noise2; single_trial_data.target{i}.neural_data.smoothed_fr(:,:,trial_nbrs(2))];
    end

    % calculate difference in firing rate between trials
    sfr_noise_diff  = (sfr_noise - sfr_noise2)/sqrt(2*min(nbr_trials_p_tgt));

    % and add a time vector for compatibility with the dim reduction code
    aux_sfr_noise_diff  = [single_trial_data.target{end}.t', sfr_noise_diff];


    % Do PCA of the noise
    dim_red_noise{j} = dim_reduction( aux_sfr_noise_diff, 'pca' );
end

% get the amount of variance explained by the first n noise components
eigenv_noise        = cell2mat(cellfun(@(x) x.eigen, dim_red_noise,'UniformOutput',false));
scree_noise         = zeros(size(eigenv_noise,1),nbr_iter);
for j = 1:nbr_iter
    scree_noise(:,j) = cumsum(eigenv_noise(:,j))/sum(eigenv_noise(:,j));
end

% and turn it into a histogram
hist_x              = 0:0.01:1;
hist_scree_noise    = zeros(length(hist_x)-1,size(eigenv_noise,1));
for j = 1:size(eigenv_noise,1)
    hist_scree_noise(:,j) = histcounts(scree_noise(j,:),hist_x)/nbr_iter;
end

% find mean noise variances
mean_noise_var      = mean(cumsum(eigenv_noise),2)/mean(sum(eigenv_noise),2);
mean_noise_eign     = mean(eigenv_noise,2);

% find 99 % limit (as in Lalazar et al., PLoC Comp Biol, 2016)
noise_var_99        = zeros(1,nbr_chs);
for i = 1:nbr_chs
    noise_var_99(i) = hist_x(find(cumsum(hist_scree_noise(:,i))>0.99,1));
end
noise_eigen_99      = prctile(eigenv_noise,99,2);



% The cutoff is when the cumulative variance explained by the real
% eigenvalues is greater than what can be proven to be not noise 
scree_trial_avg = cumsum(dim_red_trial_rel_avg.eigen)/sum(dim_red_trial_rel_avg.eigen);
% dims = find(scree_trial_avg./(alpha*(1-cumsum(noise_eigen_99))) > 1,1,'first');
dims = find( scree_trial_avg ./ (alpha*(1-cumsum(noise_eigen_99))) > 1,1,'first');


% -------------------------------------------------------------------------
% Return variables

noise_floor.dims = dims;
noise_floor.scree_trial_avg = scree_trial_avg;
noise_floor.noise_eigen_99 = noise_eigen_99;



% -------------------------------------------------------------------------
% plot both scree plots to compare

% % plot noise distribution
% figure,imagesc(1:nbr_chs,hist_x(1:end-1),hist_scree_noise),colormap('hot')
% set(gca,'TickDir','out'),set(gca,'FontSize',16);
% xlabel('Number components'), ylabel('Explained noise var.')
% 
% 
% 
% % plot scree plots neural data (different preprocessings) and noise
% figure
% hold on
% plot(cumsum(dim_red_trial_rel.eigen)/sum(dim_red_trial_rel.eigen),'linewidth',3,'marker','d','color','k')
% plot(cumsum(dim_red_trial_rel_avg.eigen)/sum(dim_red_trial_rel_avg.eigen),'linewidth',3,'marker','d','color','c')
% plot(mean_noise_var,'linewidth',3,'marker','d','color','r')
% plot(noise_var_99,'linewidth',3,'marker','d','color',[1 .6 0])
% legend('all individual trials','trial-averaged','mean noise var.','99-perc. noise var.',...
%     'Location','SouthEast'); legend boxoff;
% xlim([0 length(dim_red_trial_rel_avg.eigen)+1]),ylim([0 1])
% xlabel('Number components'), ylabel('Explained var.')
% set(gca,'TickDir','out'),set(gca,'FontSize',16);
% 
% 
% % plot eigenvalue distribution
% figure
% hold on
% plot(dim_red_trial_rel.eigen,'linewidth',3,'marker','d','color','k')
% plot(dim_red_trial_rel_avg.eigen,'linewidth',3,'marker','d','color','c')
% plot(mean_noise_eign,'linewidth',3,'marker','d','color','r')
% plot(noise_eigen_99,'linewidth',3,'marker','d','color',[1 .6 0])
% legend('all individual trials','trial-averaged','mean noise var.','99-perc. noise var.',...
%     'Location','NorthEast'); legend boxoff;
% xlim([0 length(dim_red_trial_rel_avg.eigen)+1])
% xlabel('Number components'), ylabel('variance')
% set(gca,'TickDir','out'),set(gca,'FontSize',16);
% 
