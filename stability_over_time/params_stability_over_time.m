%
% Stability over time parameters
%


% -------------------------------------------------------------------------
% Data preprocessing

% "Unsort" the sorted neurons?
pars.unsorted_yn                = true; 
% Only use electrodes that have units in all the sessions (aonly for
% unsorted_yn == true)
% pars.only_common_elecs          = false; % NO LONGER AVAILABLE as an option; there are very few channels that are stable across all sessions

% Gaussian kernel for smoothing
pars.kernel_SD                  = 0.05;

% "Downsampling rate": nbr of bins that will be combined
pars.n_bins_downs               = 3;

% Window start & end --if idx_end is empty: the duration of the shortest trial 
% WATCH OUT -- THIS IS AFTER DOWNSAMPLING SO I NEED TO UPDATE IT BASED ON
% N_BINS_DOWNS
switch lower(pars.monkey)
    case {'chips'} % S1
        pars.mani_dims = 1:8; % Neural modes to use
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 19};
    case {'han'} % S1
        pars.mani_dims = 1:8;  % Neural modes to use
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 22};
    case {'chewie','chewie2','mihili','mrt','jaco'}
        switch pars.spiking_inputs{1}
            case 'M1_spikes'
                pars.mani_dims = 1:10; % Neural modes to use
%                 pars.idx_start  = {'idx_movement_on',-2}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
%                 pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
                pars.idx_start  = {'idx_movement_on',-4}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
                pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
            case 'PMd_spikes'
                pars.mani_dims = 1:16; % Neural modes to use
                pars.idx_start  = {'idx_movement_on',-15};
                pars.idx_end    = {'idx_movement_on',0};
        end
end

% -------------------------------------------------------------------------
% To remove bad units
pars.remove_bad_neurons                 = true; % For sorted neurons it always does it
% Minimum FR per channel
pars.bad_neuron_params.min_fr           = 0.1;

% if strcmpi(pars.monkey,'chewie2'), pars.bad_neuron_params.min_fr = 0.5; end

% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;

% -------------------------------------------------------------------------
% To remove bad trials --monkey was distracted and RT was too long

% % done when data still are in 10 ms bins
% pars.bad_trial_params_PMd               = struct( ...
%                                         'ranges',{{ ...
%                                         'idx_go_cue','idx_movement_on',[0 60]; ...
%                                         'idx_target_on','idx_go_cue',[65 Inf]}});

% % done when data still are in 10 ms bins
% % SEE IF WE NEED TO CHECK THAT THE MONKEY HAD SOME REACTION TIME !!!
% pars.bad_trial_params_M1                = struct( ...
%                                         'ranges',{{ ...
%                                         'idx_movement_on','idx_trial_end',[42 Inf]}});

% % done when data still are in 10 ms bins
% pars.bad_trial_params_S1                = struct( ...
%                                         'ranges',{{ ...
%                                         'idx_movement_on','idx_trial_end',[36 Inf]}});                                    
                                    
% -------------------------------------------------------------------------
% Parameters to align latent activity
pars.align_latent_params.xval_yn        = false;
pars.align_latent_params.n_folds        = 6;

pars.align_latent_params.method         = 'cca'; % 'cca' 'procrustes'

pars.align_latent_params.n_shuff        = 100; % n shuffles for within day ceiling

pars.align_latent_params.signals        = [pars.spiking_inputs{1}(1:find(pars.spiking_inputs{1}=='_')) 'pca'];
pars.align_latent_params.mani_dims      = pars.mani_dims;

pars.align_latent_params.save_fig       = true;

pars.align_latent_params.top_lv_plot    = ceil(length(pars.mani_dims)/3);

% -------------------------------------------------------------------------
% Decoder settings -to predict behavior. In the current version, only done
% for M1 and S1, since for PMd we classify target location basde on
% preparatory activity

% Inputs and Outputs
pars.decoder_params.out                 = 'vel';
pars.decoder_params.in                  = 'aligned_data'; % 'aligned_data'; 'unaligned_data'; 'spikes' (it always does spikes afterwards)

% Bins for decoder history
pars.decoder_params.hist_bins           = 3;

% Lag of the kinematics we want to predict, only for S1 (because we predict
% past kinematics)
pars.decoder_params.lag_kin_S1          = 0; % 0.03; % (s)

% Folds for multi-fold cross-validation
pars.decoder_params.n_folds             = 6; 

% A couple other slightly redundant definitions
pars.decoder_params.manifold            = [pars.spiking_inputs{1}(1:end-7) '_pca'];
pars.decoder_params.mani_dims           = pars.mani_dims;

pars.decoder_params.unsort_chs_to_pred = true; % "unsort" channels to predict

pars.decoder_params.save_fig            = true;

% -------------------------------------------------------------------------
% Classifier settings -to predict target location based on the preparatory
% activity

% Method
pars.class_params.method                = 'Bayes'; % 'NN' 'Bayes'

% Inputs and Outputs
pars.class_params.out                   = 'target_direction';
% will do spikes after this one, by default
pars.class_params.in                    = 'aligned_data'; % 'aligned_data'; 'unaligned_data'; 'spikes'

% History?
pars.class_params.hist_bins             = 0;

% Folds for multi-fold cross-validation
pars.class_params.n_folds               = 100;
pars.class_params.num_test_trials       = 1; % number of trials per target for test set

% A couple other slightly redundant definitions
pars.class_params.manifold              = [pars.spiking_inputs{1}(1:end-7) '_pca'];
pars.class_params.mani_dims             = pars.mani_dims;

pars.class_params.idx_start             = pars.idx_start;
pars.class_params.idx_end               = pars.idx_end;
pars.class_params.idx_start_classify    = {'idx_go_cue',-13};
pars.class_params.idx_end_classify      = {'idx_go_cue',2};

pars.class_params.save_fig              = false;   
    
    
    
% -------------------------------------------------------------------------
% To estimate the "stability" of the behavior
pars.stab_behav.signal                  = 'vel';
pars.stab_behav.trial_avg               = false;
pars.stab_behav.save_fig                = true;
