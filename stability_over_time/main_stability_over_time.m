
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability of latent activity over time --main function

clear, close all



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OPTIONS
%


% -------------------------------------------------------------------------
% What data to use

pars.monkey         = 'chewie2'; % 'chewie2'; 'chewie'; 'mihili'; 'han'; 'chips'
pars.spiking_inputs = {'M1_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}; {'S1_spikes'}

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

% Neural modes to use
pars.mani_dims                  = 1:8; % 1:10; 'estimate';

% "Downsampling rate": nbr of bins that will be combined
pars.n_bins_downs               = 3;

% Window start & end --if idx_end is empty: the duration of the shortest trial 
% WATCH OUT -- THIS IS AFTER DOWNSAMPLING SO I NEED TO UPDATE IT BASED ON
% N_BINS_DOWNS
switch pars.monkey
    case {'chips'} % S1
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 19};
    case {'han'} % S1
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 22};
    case {'chewie','chewie2','mihili','mrt'}
        switch pars.spiking_inputs{1}
            case 'M1_spikes'
%                 pars.idx_start  = {'idx_movement_on',-2}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
%                 pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
                pars.idx_start  = {'idx_movement_on',-4}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
                pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
            case 'PMd_spikes'
                pars.idx_start  = {'idx_go_cue',-13}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
                pars.idx_end    = {'idx_go_cue',2}; % {''}; % {'idx_go_cue',18}
%                 pars.idx_start  = {'idx_movement_on',-15}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
%                 pars.idx_end    = {'idx_movement_on',0}; % {''}; % {'idx_go_cue',18}
        end
end
        


% -------------------------------------------------------------------------
% To remove bad units
pars.remove_bad_neurons                 = true; % For sorted neurons it always does it
% Minimum FR per channel
pars.bad_neuron_params.min_fr           = 0.1;

% if strcmp(pars.monkey,'chewie2'), pars.bad_neuron_params.min_fr = 0.5; end

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
pars.class_params.in                    = 'aligned_data'; % 'aligned_data'; 'unaligned_data'; 'spikes' (it always does spikes afterwards)

% Bins for decoder history
pars.class_params.win_size              = 0.300; % if empty, the entire window
% History?
pars.class_params.hist_bins             = 0;

% Z-score?
pars.class_params.zsc                   = true;

% Folds for multi-fold cross-validation
pars.class_params.n_folds               = 5; 

% A couple other slightly redundant definitions
pars.class_params.manifold              = [pars.spiking_inputs{1}(1:end-7) '_pca'];
pars.class_params.mani_dims             = pars.mani_dims;

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


% -------------------------------------------------------------------------
% Datafiles for each monkey


files_chewie        = { ...
                        'Chewie_CO_VR_2016-09-09.mat', ...
                        'Chewie_CO_VR_2016-09-12.mat', ...
                        'Chewie_CO_VR_2016-09-14.mat', ...
                        'Chewie_CO_FF_2016-09-15.mat', ...
                        'Chewie_CO_FF_2016-09-19.mat', ...
                        'Chewie_CO_FF_2016-09-21.mat', ...
                        'Chewie_CO_FF_2016-09-23.mat', ...
                        'Chewie_CO_VR_2016-09-29.mat', ...
                        'Chewie_CO_FF_2016-10-05.mat', ...
                        'Chewie_CO_VR_2016-10-06.mat', ...
                        'Chewie_CO_FF_2016-10-07.mat', ...
                        'Chewie_CO_FF_2016-10-11.mat', ...
                        'Chewie_CO_FF_2016-10-13.mat', ...
                        'Chewie_CO_CS_2016-10-14.mat', ...
                        'Chewie_CO_CS_2016-10-21.mat' ...
                        };
 % M1: Get rid of a session with much longer RT than the others, and with
 % a bunch of trials with idx_movement_onset = NaN
if strcmp(pars.spiking_inputs{1},'M1_spikes') && strcmp(pars.monkey,'chewie')
    files_chewie( strncmpi(files_chewie,'Chewie_CO_FF_2016-09-23.mat',...
        length('Chewie_CO_FF_2016-09-23.mat')) ) = [];
end            
    

files_chewie2       = { ...
                        'Chewie_CO_VR_2013-10-03.mat', ...
                        'Chewie_CO_FF_2013-10-22.mat', ...
                        'Chewie_CO_FF_2013-10-23.mat', ...
                        'Chewie_CO_FF_2013-10-31.mat', ...
                        'Chewie_CO_FF_2013-11-01.mat', ...
                        'Chewie_CO_FF_2013-12-03.mat', ...
                        'Chewie_CO_FF_2013-12-04.mat', ...
                        'Chewie_CO_VR_2013-12-19.mat', ...
                        'Chewie_CO_VR_2013-12-20.mat', ...
                        'Chewie_CO_CS_2015-03-09.mat', ...
                        'Chewie_CO_CS_2015-03-11.mat', ...
                        'Chewie_CO_CS_2015-03-12.mat', ...
                        'Chewie_CO_CS_2015-03-13.mat', ...
                        'Chewie_CO_CS_2015-03-19.mat', ...
                        'Chewie_CO_FF_2015-06-29.mat', ...
                        'Chewie_CO_FF_2015-06-30.mat', ...
                        'Chewie_CO_FF_2015-07-01.mat', ...
                        'Chewie_CO_FF_2015-07-03.mat', ... 
                        'Chewie_CO_FF_2015-07-06.mat', ...
                        'Chewie_CO_FF_2015-07-07.mat', ...
                        'Chewie_CO_FF_2015-07-08.mat', ...
                        'Chewie_CO_VR_2015-07-09.mat', ...
                        'Chewie_CO_VR_2015-07-10.mat', ...
                        'Chewie_CO_VR_2015-07-13.mat', ...
                        'Chewie_CO_VR_2015-07-14.mat', ...
                        'Chewie_CO_VR_2015-07-15.mat', ...
                        'Chewie_CO_VR_2015-07-16.mat', ...
                        'Chewie_CO_CS_2015-11-03.mat', ...
                        'Chewie_CO_CS_2015-11-04.mat', ...
                        'Chewie_CO_CS_2015-11-06.mat', ...
                        'Chewie_CO_GR_2015-11-09.mat', ...
                        'Chewie_CO_GR_2015-11-10.mat', ...
                        'Chewie_CO_GR_2015-11-12.mat', ...
                        'Chewie_CO_GR_2015-11-13.mat', ...
                        'Chewie_CO_GR_2015-11-16.mat', ...
                        'Chewie_CO_GR_2015-11-17.mat', ...
                        'Chewie_CO_VR_2015-11-19.mat', ...
                        'Chewie_CO_GR_2015-11-20.mat', ...
                        'Chewie_CO_VR_2015-12-01.mat', ...
                        'Chewie_CO_VR_2015-12-03.mat', ...
                        'Chewie_CO_VR_2015-12-04.mat' ...
                        };

% Delete some sessions without an instructed delya
if strcmp(pars.spiking_inputs{1},'M1_spikes') && strcmp(pars.monkey,'chewie2')
    files_chewie2( strncmpi(files_chewie2,'Chewie_CO_VR_2013-10-03.mat',...
        length('Chewie_CO_VR_2013-10-03.mat')) ) = [];
    files_chewie2( strncmpi(files_chewie2,'Chewie_CO_FF_2013-10-22.mat',...
        length('Chewie_CO_FF_2013-10-22.mat')) ) = [];
    files_chewie2( strncmpi(files_chewie2,'Chewie_CO_FF_2013-10-23.mat',...
        length('Chewie_CO_FF_2013-10-23.mat')) ) = [];
    files_chewie2( strncmpi(files_chewie2,'Chewie_CO_FF_2013-10-31.mat',...
        length('Chewie_CO_FF_2013-10-31.mat')) ) = [];
    files_chewie2( strncmpi(files_chewie2,'Chewie_CO_FF_2013-11-01.mat',...
        length('Chewie_CO_FF_2013-11-01.mat')) ) = [];
end            
                    

files_mihili        = { ...
                        'Mihili_CO_FF_2014-02-03.mat', ...
                        'Mihili_CO_FF_2014-02-17.mat', ...
                        'Mihili_CO_FF_2014-02-18.mat', ...
                        'Mihili_CO_VR_2014-03-03.mat', ...
                        'Mihili_CO_VR_2014-03-04.mat', ...
                        'Mihili_CO_VR_2014-03-06.mat', ...
                        'Mihili_CO_FF_2014-03-07.mat', ...
                        'Mihili_CO_CS_2014-09-29.mat', ...
                        'Mihili_CO_CS_2014-12-03.mat', ...
                        'Mihili_CO_CS_2015-05-11.mat', ...
                        'Mihili_CO_FF_2015-06-10.mat', ...
                        'Mihili_CO_FF_2015-06-11.mat', ...
                        'Mihili_CO_FF_2015-06-15.mat', ...
                        'Mihili_CO_FF_2015-06-16.mat', ...
                        'Mihili_CO_FF_2015-06-17.mat', ...
                        'Mihili_CO_VR_2015-06-23.mat', ...
                        'Mihili_CO_VR_2015-06-25.mat', ...
                        'Mihili_CO_VR_2015-06-26.mat' ...
                        };

% M1: Get rid of a session with very few channels with neurons (9!)
if strcmp(pars.spiking_inputs{1},'M1_spikes') && strcmp(pars.monkey,'mihili')
    files_mihili( strncmpi(files_mihili,'Mihili_CO_CS_2014-12-03.mat',...
        length('Mihili_CO_CS_2014-12-03.mat')) ) = [];
end
% PMd: Get rid of two sessions without PMd data
if strcmp(pars.spiking_inputs{1},'PMd_spikes') && strcmp(pars.monkey,'mihili')
    files_mihili( strncmpi(files_mihili,'Mihili_CO_CS_2014-06-26.mat',...
        length('Mihili_CO_CS_2014-06-26.mat')) ) = [];
    files_mihili( strncmpi(files_mihili,'Mihili_CO_CS_2014-06-27.mat',...
        length('Mihili_CO_CS_2014-06-27.mat')) ) = [];
end
                    
files_mrt           = { ...
                        'MrT_CO_FF_2013-08-19.mat', ...
                        'MrT_CO_FF_2013-08-21.mat', ...
                        'MrT_CO_FF_2013-08-23.mat', ...
                        'MrT_CO_VR_2013-09-03.mat', ...
                        'MrT_CO_VR_2013-09-05.mat', ...
                        'MrT_CO_VR_2013-09-09.mat' ...
                        };
            
files_han           = { 'Han_COactpas_2017-10-24.mat', ...
                        'Han_COactpas_2017-10-30.mat', ...
                        'Han_COactpas_2017-10-31.mat', ...
                        'Han_COactpas_2017-11-03.mat', ...
                        'Han_COactpas_2017-11-16.mat', ...
                        'Han_COactpas_2017-11-20.mat', ...
                        'Han_COactpas_2017-11-21.mat', ...
                        'Han_COactpas_2017-11-22.mat', ...
                        'Han_COactpas_2017-11-27.mat', ...
                        'Han_COactpas_2017-11-28.mat', ...
                        'Han_COactpas_2017-11-29.mat', ...
                        'Han_COactpas_2017-12-01.mat', ...
                        'Han_COactpas_2017-12-04.mat', ...
                        'Han_COactpas_2017-12-07.mat' ...
                        }; 

files_chips         = { ...
                        'Chips_20151113_TD_nosort_notrack_noemg.mat', ...
                        'Chips_20151117_TD_nosort_notrack_noemg.mat', ...
                        'Chips_20151120_TD_nosort_notrack_noemg.mat', ...
                        'Chips_20151201_TD_nosort_notrack_noemg.mat', ...
                        'Chips_20151204_TD_nosort_notrack_noemg.mat', ...
                        'Chips_20151211_TD_nosort_notrack_noemg.mat' ...
                        };


% -------------------------------------------------------------------------
% Go to the path where the data are, and update the file list if necessary
% (i.e., if using "unsorted" data)

here                = pwd;

% If we want to "unsort", choose the appropriate file names
switch pars.monkey
    % Chewie M1-PMd
    case 'chewie'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        files = files_chewie;
    % Chewie M1 only implant
    case 'chewie2'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        files = files_chewie2;
    % Mihili M1-PMd
    case 'mihili'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili')
        files = files_mihili;
    % MrT PMd
    case 'mrt'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/MrT')
        files = files_mrt;
    % Han S1
    case 'han'
        if ~pars.unsorted_yn
            warning('Data for Han are not sorted'); 
            warning('Loading threshold crossings'); 
            pause(3);
            pars.unsorted_yn = true;
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
        end
        files = files_han;
    % Chips S1
    case 'chips'
        if ~pars.unsorted_yn
            warning('Data for Chips are not sorted'); 
            warning('Loading threshold crossings'); 
            pause(3);
            pars.unsorted_yn = true;
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chips')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chips')
        end
        files = files_chips; 
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
    if sum(strcmp(pars.monkey,{'mihili','chewie','chewie2','mrt'}))


        if strcmp(pars.spiking_inputs{1},'PMd_spikes')
        
            [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@stripSpikeSorting},...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)} , ...
                            {@trimTD,pars.idx_start,pars.idx_end});
                        
        elseif strcmp(pars.spiking_inputs{1},'M1_spikes')
            
            if ~strcmp(pars.monkey,'chewie2')
                
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
            elseif strcmp(pars.monkey,'chewie2')
         
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
    elseif sum(strcmp(pars.monkey,{'chips','han'}))
        
        [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'result','R'}, ...
                            {@stripSpikeSorting});
    end
end

% Keep PCA pars, get rid of the rest of the extra outputs
if strcmp(pars.monkey,'chewie2')
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end};
    end
elseif ~( strcmp(pars.monkey,'han') || strcmp(pars.monkey,'chips') )
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
if strcmp(pars.monkey,'han')
    
    % ---------------------------------------------------------------------
    % there are trials with targets = NaN -> exclude them
    master_td       = master_td(~isnan([master_td.target_direction]));
    
    % ---------------------------------------------------------------------
    % in some sessions he had 4 targets and in others. Only keep the common
    % targets
    all_targets     = unique([master_td.target_direction]);
    targets_session = nan(n_sessions,length(all_targets));
    for s = 1:n_sessions
        this_s = getTDidx(master_td,{'date',meta.sessions{s}});
        u_targets = unique([master_td(this_s).target_direction]);
        targets_session(s,1:length(u_targets)) = u_targets;
    end
    
    is_target = zeros(n_sessions,length(all_targets));
    for s = 1:n_sessions
        is_target(s,:) = ismember(all_targets,targets_session(s,:));
    end
    
    % Targets that were present in all sessions
    targets_always_present = all_targets(sum(is_target,1)==n_sessions);
    
    % Only keep trials to targets that were present in all the sessions
    [~,master_td] = getTDidx(master_td,{'target_direction',targets_always_present});

    
    clear unique_targets_session is_target targets_always_present u_targets targets_session this_s all_targets;
% -- And for Chips, to exclude NaN targets
elseif strcmp(pars.monkey,'chips')
     % there are trials with targets = NaN -> exclude them
    master_td       = master_td(~isnan([master_td.target_direction]));
end


meta.targets        = unique([master_td.target_direction]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DO PCA FOR THE S1 MONKEYS --because we needed the fixes above
%


if strcmp(pars.monkey,'han') || strcmp(pars.monkey,'chips')
   
    % ---------------------------------------------------------------------
    % Calculate movement onset
    master_td       = getMoveOnsetAndPeak(master_td,struct('start_idx','idx_go_cue','end_idx','idx_trial_end'));
  
%     % ---------------------------------------------------------------------
%     % Remove bad trials
% %     % I CAN'T MAKE THIS KLJSAD WORK WELL
% %     master_td       = removeBadTrials( master_td, struct('range',pars.bad_trial_params_S1,'remove_nan_idx',false) );
%     RT              = [master_td.idx_movement_on] - [master_td.idx_go_cue];
%     MT              = [master_td.idx_trial_end] - [master_td.idx_movement_on];
%     master_td( MT<= ( (pars.idx_end{2}-pars.idx_start{2}+1) * pars.n_bins_downs ) ) = [];
    
    
    % ---------------------------------------------------------------------
    % Downsample
    master_td       = binTD( master_td, pars.n_bins_downs );
    
    % ---------------------------------------------------------------------
    % Remove bad neurons, square root transform and smooth for PCA
    
    new_master_td   = [];
    for s = 1:n_sessions
        
        this_s      = getTDidx( master_td, {'date',meta.sessions{s}} );
        
        % remove bad units
        this_td     = removeBadNeurons( master_td(this_s), pars.bad_neuron_params );
        
        % square root transform
        this_td     = sqrtTransform( this_td, pars.spiking_inputs );
        
        % smooth
        this_td     = smoothSignals( this_td, struct(...
                            'signals',pars.spiking_inputs, ...
                            'calc_fr', true, ...
                            'kernel_SD', pars.kernel_SD) );
                        
        new_master_td = [new_master_td, this_td];
    end        
    
    master_td       = new_master_td; 
    clear new_master_td;
    
    % ---------------------------------------------------------------------
    % trim the trials to idx_target_on : idx_trial_end, the window for PCA
    master_td       = trimTD( master_td, {'idx_target_on',0}, {'idx_trial_end',0} );
        
    % ---------------------------------------------------------------------
    % Do PCA
    new_td          = [];
    for s = 1:n_sessions
        
        this_s = getTDidx(master_td,{'date',meta.sessions{s}});
        
        [this_td, pca_info(s)] = getPCA( master_td(this_s), struct('signals',pars.spiking_inputs) );
        
        new_td = [new_td, this_td];
    end
    
    master_td       = new_td;
        
    % ---------------------------------------------------------------------
    % TRIM THE TRIALS to the ANALYSIS WINDOW
    master_td           = trimTD( master_td, pars.idx_start, pars.idx_end );
    
    % GET RID OF TARGETS THAT ARE TOO SHORT BRUTE FORCE BECAUSE I HATE THIS
    too_short           = arrayfun(@(x) size(x.pos,1), master_td) < (pars.idx_end{2} - pars.idx_start{2} + 1);
    disp([num2str(sum(too_short)) ' Trials are outside the window']);
    master_td(too_short) = [];
end


% -------------------------------------------------------------------------
% HACK FOR CHEWIE M1 ONLY DATA: in a few trials, the time between movement
% onset and reward is very short --discard those

% if strcmp(pars.monkey,'chewie2')
%     
%     % Analysis window duration
%     analysis_win = pars.idx_end{2} - pars.idx_start{2};
%     
%     % Get window duration for each trial, assuming it goes from movement
%     % onset until trial end since this is M1
%     if strcmp(pars.idx_start{1},'idx_movement_on') || strcmp(pars.idx_end{1},'idx_movement_on')
%         
%         win_p_trial = arrayfun( @(x) x.idx_trial_end-x.idx_movement_on, master_td);
%         trials_discard = win_p_trial <= analysis_win;
%         
%         master_td(trials_discard) = [];
%         
%         disp(['Discarded ' num2str(sum(trials_discard)) ' because of their very short MT']);
%     else
%         error('This is for M1! Fix the analysis window so it is expressed as fcn of movement onset, or add check to this code!');
%     end
% end


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
% ESTIMATE DIMENSIONALITY
%

% % NEEDTOFIX: SINCE THE BINS ARE PRETTY WIDE AND WE HAVE ONLY 4 TRIALS AND A
% % LOT OF NEURONS, THE PCA FUNCTION DOESN'T WORK 

% if strcmp(pars.spiking_inputs{1},'M1_spikes') % || strcmp(pars.spiking_inputs{1},'PMd_spikes')
%  
%     for s = 1:n_sessions
% 
%         this_s          = getTDidx(master_td,{'date',meta.sessions{s}});
%         dims(s)         = estimateDimensionality(master_td(this_s),struct(...
%                                 'signals',pars.spiking_inputs,...
%                                 'alpha',0.95,...
%                                 'num_iter',1000));
%     end
%     figure, histogram(dims,0:10,'facecolor','k')
%     set(gca,'TickDir','out','FontSize',14), box off
%     xlabel('Dimensionality'),ylabel('Counts'),
%     title([pars.monkey ' - ' pars.spiking_inputs{1}(1:end-7)])
% end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COMPARE KINEMATICS
%

if ~strcmp(pars.spiking_inputs{1},'PMd_spikes')
    corr_kin        = comp_behavior( master_td, pars.stab_behav );
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ALIGN LATENT ACTIVITY OVER SESSIONS
%

align_results       = align_latent_activity( master_td, pars.align_latent_params );

% % get within-day ceiling
% if strcmp(pars.monkey,'chewie2')
%     master_td_all_trials = master_td;
% end
within_day_align_results = align_latent_activity_within_day( master_td_all_trials, pars.align_latent_params );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DECODE KINEMATICS TRAINING BUILDING A DECODER ON ONE DAY AND TESTING IT
% ON ANOTHER USING THE ALIGNED LATENT ACTIVITY
%
% only do decoding for S1 / M1 for PMd we'll classify target direction

if strcmp(pars.spiking_inputs{1},'M1_spikes') || strcmp(pars.spiking_inputs{1},'S1_spikes')
    
    dec_results     = decode_across_days( master_td, pars.decoder_params );
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DECODE KINEMATICS TRAINING BUILDING A DECODER ON ONE DAY AND TESTING IT
% ON ANOTHER USING THE ALIGNED LATENT ACTIVITY
%
% only do decoding for S1 / M1 for PMd we'll classify target direction

if strcmp(pars.spiking_inputs{1},'M1_spikes') || strcmp(pars.spiking_inputs{1},'S1_spikes')
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

% if strcmp(pars.spiking_inputs{1},'PMd_spikes')
%    
%     % Do for many window sizes
%     for w = 1:16
%         
%         pars.class_params.win_size = master_td(1).bin_size * w;
%         
%         clas_results    = classify_across_days( master_td, pars.class_params );
%     end
% end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PLOTTING



% Plot similarity behavior
if strcmp(pars.spiking_inputs{1},'M1_spikes') || strcmp(pars.spiking_inputs{1},'S1_spikes')
    SOT_Fig_stability_behavior( corr_kin, pars );
end


% Plot aligned latent activity and similarity over days
SOT_Fig_3_aligned_latent_activity( master_td, align_results, meta, pars.align_latent_params, within_day_align_results );


% Plot decoding results
if strcmp(pars.spiking_inputs{1},'M1_spikes') || strcmp(pars.spiking_inputs{1},'S1_spikes')
    
    SOT_Fig_decoding( dec_results, dec_spike_results, pars );
elseif strcmp(pars.spiking_inputs{1},'PMd_spikes')
    
    
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

