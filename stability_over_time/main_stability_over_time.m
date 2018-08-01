
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability of latent activity over time --main function

clear, close all

% get the directories for the computer running the code
[pars.save_dir, pars.data_dir] = get_computer_paths();


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OPTIONS
%


% -------------------------------------------------------------------------
% What data to use

pars.monkey         = 'Chips'; % 'chewie2'; 'chewie'; 'mihili'; 'han'; 'chips'; 'jaco'
pars.spiking_inputs = {'S1_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}; {'S1_spikes'}

% Sesssions to discard if any
pars.sessions_discard = []; %6:14; %[12 13 14];

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
            case 'M1s_spikes'
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOAD THE DATA
%
monkey_file_lists;

% -------------------------------------------------------------------------
% Go to the path where the data are, and update the file list if necessary
% (i.e., if using "unsorted" data)

here                = pwd;

% If we want to "unsort", choose the appropriate file names
switch lower(pars.monkey)
    % Chewie M1-PMd
    case 'chewie'
        %cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        files = files_chewie;
    % Chewie M1 only implant
    case 'chewie2'
        %cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        files = files_chewie2;
    case 'jaco'
        files = files_jaco;
    % Mihili M1-PMd
    case 'mihili'
        %cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili')
        files = files_mihili;
    % MrT PMd
    case 'mrt'
        %cd('/Users/juangallego/Documents/NeuroPlast/Data/MrT')
        files = files_mrt;
    % Han S1
    case 'han'
        if ~pars.unsorted_yn
            warning('Data for Han are not sorted'); 
            warning('Loading threshold crossings'); 
            pause(3);
            pars.unsorted_yn = true;
            %cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
        else
            %cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
        end
        files = files_han;
    % Chips S1
    case 'chips'
        if ~pars.unsorted_yn
            warning('Data for Chips are not sorted'); 
            warning('Loading threshold crossings'); 
            pause(3);
            pars.unsorted_yn = true;
            %cd('/Users/juangallego/Documents/NeuroPlast/Data/Chips')
        else
            %cd('/Users/juangallego/Documents/NeuroPlast/Data/Chips')
        end
        files = files_chips; 
end

% turn the files into a full filename
for i = 1:length(files)
    files{i} = fullfile(pars.data_dir,pars.monkey,files{i});
end


% -------------------------------------------------------------------------
% Load the data

% Sorted data
if ~pars.unsorted_yn
    
    % This loads and already preprocesses and does dimensionality reduction
    % 1) Loads Baseline trials only (no force field or visuomotor
    % adaptation); 2) Loads only successful trials; 3) Square
    % root-transforms, smooths and does PCA (on the entire trial); 4)
    % Downsamples the data 
    [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@removeBadTrials,pars.bad_trial_params}, ...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTrials,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)} );  
    
% Unsort the units                        
else      
    
    % % This loads and only does a little bit of preprocessing, with no dim
    % reduction: 1) Loads Baseline trials only (no force field or visuomotor
    % adaptation); 2) Loads only successful trials; 3) Merges all units in
    % the same channel for Mihili and Chewie, but not for Han
    if sum(strcmpi(pars.monkey,{'mihili','chewie','chewie2','mrt'}))


        if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
        
            [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@stripSpikeSorting},... % @removeBadTrials, ...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)});
                        
        elseif strcmpi(pars.spiking_inputs{1},'M1_spikes')
            
            if ~strcmpi(pars.monkey,'chewie2')
                
                [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@stripSpikeSorting},...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)}, ...
                            {@trimTD,pars.idx_start,pars.idx_end});
                        
            % For Chewie2: we have to remove some trials for which end time
            % goes beyong trial time, and some others for which requested
            % start time is <1
            elseif strcmpi(pars.monkey,'chewie2')
         
                % request minimum movement time and delay between trial
                % start and movement onset
                min_MT = (pars.idx_end{2}+1)*pars.n_bins_downs;
                min_WT = 20; % (abs(pars.idx_start{2})+2)*pars.n_bins_downs;
                
                pars.bad_trial_params = struct( 'ranges',{{ ...
                                        'idx_movement_on','idx_trial_end',[min_MT Inf]; ...
                                        'idx_target_on','idx_movement_on',[min_WT Inf]}});
                
                [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@stripSpikeSorting},...
                            {@removeBadTrials,pars.bad_trial_params}, ...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)}); 
                
%                 % I have to remove params again, because I'm dumb and can't
%                 % do it in one step above -- but it works
%                 bp = struct('ranges',{{'idx_target_on','idx_movement_on',[5 Inf]}},'remove_nan_idx',false);
%                 master_td = removeBadTrials(master_td,bp);
                
                % Now trim
                master_td = trimTD( master_td, pars.idx_start, pars.idx_end );
            end 
        end

    % FOR HAN and CHIPS
    elseif sum(strcmpi(pars.monkey,{'chips','han'}))
        
        [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'result','R'}, ...
                            {@stripSpikeSorting});
    end
end

% Keep PCA pars, get rid of the rest of the extra outputs
if strcmpi(pars.monkey,'chewie2')
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end};
    end
elseif strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end};
    end
elseif ~( strcmpi(pars.monkey,'han') || strcmpi(pars.monkey,'chips') )
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end-1};
    end
clear pars_td s;
end


% go back to where you were path-wise
cd(here);
clear files* here;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET METADATA
%


% -------------------------------------------------------------------------
% GET THE SESSIONS 
meta.sessions       = unique({master_td.date});
n_sessions          = length(meta.sessions);


% Sometimes the sessions are not sorted by time --fix that here by
% resorting the trials in master_td
[~, i_sort]         = sort(datenum([meta.sessions]));
if sum( i_sort - 1:length(i_sort) ) > 0
   
    sorted_dates = sort( cell2mat( cellfun(@(x) datenum(x), meta.sessions, 'uni', 0) ) );
    for s = 1:n_sessions
        meta.sessions{s} = datestr(sorted_dates(s),'mm-dd-yyyy');
    end
end


% -------------------------------------------------------------------------
% DISCARD SESSIONS YOU WANT TO EXCLUDE

if ~isempty(pars.sessions_discard)
    
    for s = 1:length(pars.sessions_discard)
        this_s = getTDidx(master_td,{'date',meta.sessions{pars.sessions_discard(s)}});
        master_td(this_s) = [];
    end
    
    % Update the sessions
    meta.sessions   = unique({master_td.date});
    n_sessions      = length(meta.sessions);
    
    clear this_s;
end


% -------------------------------------------------------------------------
% GET THE TARGETS

% --For Han and Chips, we need a couple of tricks, because of errors in
% data storage
if strcmpi(pars.monkey,'han') || strcmpi(pars.monkey,'chips')

    prep_to_get_targets_S1;
end

meta.targets        = unique([master_td.target_direction]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DO PCA FOR THE S1 MONKEYS --because we needed the fixes above
%


if strcmpi(pars.monkey,'han') || strcmpi(pars.monkey,'chips')
   
    prep_and_do_PCA_S1;
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EQUALIZE THE NUMBER OF TRIALS ACROSS TARGETS AND SESSIONS
%

% save master_td to align within days, instead of using the equalized
% version
master_td_all_trials = master_td;


% Equalize
master_td           = equalNbrTrialsSessions(master_td);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COMPARE KINEMATICS
%

if ~strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    corr_kin        = comp_behavior( master_td, pars.stab_behav );
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ALIGN LATENT ACTIVITY OVER SESSIONS
%
if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    % need to trim it first here
    align_results       = align_latent_activity( trimTD(master_td, pars.idx_start, pars.idx_end), pars.align_latent_params );
    within_day_align_results = align_latent_activity_within_day( trimTD(master_td_all_trials, pars.idx_start, pars.idx_end), ...
                                pars.align_latent_params );
else
    align_results       = align_latent_activity( master_td, pars.align_latent_params );
    within_day_align_results = align_latent_activity_within_day( master_td_all_trials, pars.align_latent_params );
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DECODE KINEMATICS TRAINING BUILDING A DECODER ON ONE DAY AND TESTING IT
% ON ANOTHER USING THE ALIGNED LATENT ACTIVITY
%
% only do decoding for S1 / M1 for PMd we'll classify target direction

if strcmpi(pars.spiking_inputs{1},'M1_spikes') || strcmpi(pars.spiking_inputs{1},'S1_spikes')
    
    dec_results     = decode_across_days( master_td, pars.decoder_params );
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DECODE KINEMATICS TRAINING BUILDING A DECODER ON ONE DAY AND TESTING IT
% ON ANOTHER USING THE ALIGNED LATENT ACTIVITY
%
% only do decoding for S1 / M1 for PMd we'll classify target direction

if strcmpi(pars.spiking_inputs{1},'M1_spikes') || strcmpi(pars.spiking_inputs{1},'S1_spikes')
    in_old          = pars.decoder_params.in;
    pars.decoder_params.in = 'spikes';

    dec_spike_results = decode_across_days( master_td, pars.decoder_params );

    pars.decoder_params.in = in_old; clear in_old;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CLASSIFY TARGET LOCATION
%

if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    
    % first do aligned or unaligned
    clas_results = classify_across_days(master_td,pars.class_params.in,pars.class_params);

    % now do spikes
    clas_spike_results = classify_across_days(master_td,'spikes',pars.class_params);
    clas_spike_results.norm_perf_spike =  clas_spike_results.perf_spike./mean(clas_results.perf_within_xval2,2)*100;
    
    % now do plotting
    %%%% TO DO
    
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PLOTTING



% Plot similarity behavior
if strcmpi(pars.spiking_inputs{1},'M1_spikes') || strcmpi(pars.spiking_inputs{1},'S1_spikes')
    SOT_Fig_stability_behavior( corr_kin, pars );
end


% Plot aligned latent activity and similarity over days
if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    % need to trim it first here
    SOT_Fig_3_aligned_latent_activity( trimTD(master_td,pars), pars.save_dir, align_results, meta, pars.align_latent_params, within_day_align_results );
else
    SOT_Fig_3_aligned_latent_activity( master_td, pars.save_dir, align_results, meta, pars.align_latent_params, within_day_align_results );
end


% Plot decoding results
if strcmpi(pars.spiking_inputs{1},'M1_spikes') || strcmpi(pars.spiking_inputs{1},'S1_spikes')
    SOT_Fig_decoding( dec_results, dec_spike_results, pars );
elseif strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    SOT_Fig_classification(clas_results, clas_spike_results, pars);
end



% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % CONTROLS
% %
% 
%
% % 1. Manifold stability over time
%
% pars.align_latent_params.surrogate_type = 'surrogate-T'; % 'surrogate-C' 'surrogate-TC'
% 
% tme_align_results = align_latent_tme_control( master_td, pars.align_latent_params );
%
% 
% % 2. Principal angles over time
% 
% mani_comp_results = comp_manifolds_across_days( master_td, pars );
%
% 
% % 3. Artificially modify single unit properties - TODO
% 
% shuff_align_results = align_latent_shuffle_unit_props( master_td, pars.align_latent_params );

