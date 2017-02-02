%
% Meta wrapper to run analysis for the manifold paper
%


pool                            = gcp;


% % load data, if not in the workspace
% if ~exist('datasets','var')
%     path            = '/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction';
%     load([path filesep 'all_manifold_datasets.mat']);
% end


% -------------------------------------------------------------------------
% % 0. these data will be already preprocessed, if not:
prep_params                     = batch_preprocess_dim_red_data_params_defaults;
prep_params.bin_size            = .02;
prep_params.save_data           = false;
prep_params.norm_trial_data     = 'min_dur'; % 'min_dur' 't_warp'
prep_params.dim_red_emg         = 'nmf'; % 'nmf'
prep_params.emg_factors         = 3;
% preprocess the data
datasets                        = batch_preprocess_dim_red_data( pwd, prep_params );

% compare dim_red_emg data
if ~strcmp(prep_params.dim_red_emg,'none') && ~prep_params.dim_red_emg_across_taks 
    compare_muscle_synergy_spaces_all(datasets,prep_params.emg_factors);
end


% -------------------------------------------------------------------------
% 1. dimensionality analysis: compare variance explained with PCA, etc
dim_results                     = batch_dimensionality_analysis( datasets );


% -------------------------------------------------------------------------
% 2. angles between manifold dimensions
angle_results                   = batch_angle_btw_manifold_dims( datasets );


% -------------------------------------------------------------------------
% 3. compare the dynamics of the 'neural projections'
% load the default params...
proj_params                     = batch_compare_manifold_projs_defaults();
proj_params.nbr_shuffles_bootstrap = 10000; % instead of the default 1000
% and run the analysis
proj_results                    = batch_compare_manifold_projs( datasets, proj_params );


% -------------------------------------------------------------------------
% 4. compute and compare task potent / null spaces
opn_params                      = batch_compare_potent_null_spaces_defaults;
opn_params.time_win             = [ 0.35 1.15; 0.4 1.2; 0.3 1.1;
                                    0 0.7; 0 0.7; 0 0.7; 
                                    0.15 0.95; 0.15 0.95 ];
opn_params.dim_neural_manifold  = 10;
opn_params.dim_emg_manifold     = 6;
opn_params.neural_to_output_delay = .05;
opn_params.output               = 'dim_red_emg'; % 'dim_red_emg' 'emg'
opn_params.trial_averaged       = true;
% run the analysis
opn_results                     = batch_compare_potent_null_spaces( datasets, opn_params );




% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
  
% % OLD CODE ---
% % NEED TO HAVE LOADED A BUNCH OF THINGS WITH RUN_DIM_ANALYSIS_2.M
% [analysis.onp.dim_raw, analysis.onp.dim_summary, analysis.onp.single_trial_data ] = call_find_output_null_potent_dims_updated( ...
%         binned_data, neural_chs, chosen_emgs, labels, smoothed_FR, false);


% % define the number of muscle synergies
% nbr_muscle_synergies = 4;
% data                = cell(1,length(datasets));
% 
% % compare muscle spaces across tasks
% for d = 1:length(datasets)
%     % retrieve PCA matrices for each task
%     W_emg           = cellfun( @(x) x.target{end}.emg_data.dim_red.w , datasets{d}.stdata, ...
%                         'UniformOutput', false );
%     % get all combs of tasks
%     comb_bdfs       = nchoosek(1:length(W_emg),2);
%     nbr_comb_bdfs   = size(comb_bdfs,1);
%     
%     % do SVD of the synergy matrix
%     [~, ~, V_emg]   = cellfun( @(x) svd(x), W_emg, 'UniformOutput', false );
% 
%     pr_angles_emg   = zeros(nbr_comb_bdfs,nbr_muscle_synergies);
%     
%     % and compare their synergy spaces 
%     for p = 1:nbr_comb_bdfs
%         pr_angles_emg(p,:) = principal_angles( V_emg{comb_bdfs(p,1)}(:,1:nbr_muscle_synergies), ...
%                                 V_emg{comb_bdfs(p,2)}(:,1:nbr_muscle_synergies) );
%     end
%     
%     % store
%     data{d}.princ_angles_emg = pr_angles_emg;
% end
% 
% 
% % plot
% colors              = distinguishable_colors(length(data));
% figure, hold on
% xlim([0 nbr_muscle_synergies+1]),ylim([0 90])
% set(gca,'Tickdir','out'),set(gca,'FontSize',16)
% xlabel('dimension'),ylabel('angle (deg)')
% for d = 1:length(data)
%     plot( rad2deg(data{d}.princ_angles_emg)', 'linewidth', 2, 'color', colors(d,:) );
% end
% 
