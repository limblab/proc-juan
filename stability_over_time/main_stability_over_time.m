
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

pars.monkey         = 'chewie'; % 'chewie'; 'mihili'; 'han'
pars.spiking_inputs = {'M1_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}

% Sesssions to discard if any
pars.sessions_discard = []; %[12 13 14];

% -------------------------------------------------------------------------
% Data preprocessing

% "Unsort" the sorted neurons?
pars.unsorted_yn    = false; 
% Only use electrodes that have units in all the sessions (aonly for
% unsorted_yn == true)
pars.only_common_elecs = false;

% Gaussian kernel for smoothing
pars.kernel_SD      = 0.1;

% Neural modes to use
pars.mani_dims      = 1:10; % 1:10; 'estimate';

% Window start & end --if idx_end is empty: the duration of the shortest trial 
switch pars.monkey
    case 'han'
        pars.idx_start      = {'idx_goCueTime',round(-2*5/3)};% {'idx_goCueTime',-2}; % {'idx_movement_on',-2}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
        pars.idx_end        = {'idx_goCueTime',round(9*5/3)}; % {'idx_goCueTime',9}; % {'idx_movement_on',13}; % {''}; % {'idx_go_cue',18}
    case {'chewie','mihili'}
        pars.idx_start      = {'idx_movement_on',-2}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
        pars.idx_end        = {'idx_movement_on',13}; % {''}; % {'idx_go_cue',18}
end
        
% "Downsampling rate": nbr of bins that will be combined
pars.n_bins_downs   = 3;


% -------------------------------------------------------------------------
% To remove bad units
pars.remove_bad_neurons                 = true; % For sorted neurons it always does it
% Minimum FR per channel
pars.bad_neuron_params.min_fr           = 0.1;
% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;


% -------------------------------------------------------------------------
% Parameters to align latent activity
pars.align_latent_params.xval_yn        = false;
pars.align_latent_params.n_folds        = 6;

pars.align_latent_params.method         = 'cca'; % 'cca' 'procrustes'

pars.align_latent_params.signals        = [pars.spiking_inputs{1}(1:find(pars.spiking_inputs{1}=='_')) 'pca'];
pars.align_latent_params.mani_dims      = pars.mani_dims;

% -------------------------------------------------------------------------
% Decoder settings -to predict behavior

% Inputs and Outputs
pars.decoder_params.out                 = 'vel';
pars.decoder_params.in                  = 'aligned_data'; % 'aligned_data'; 'raw_data'

% Bins for decoder history
pars.decoder_params.hist_bins           = 3;

% Folds for multi-fold cross-validation
pars.decoder_params.n_folds             = 6; 

% A couple other slightly redundant definitions
pars.decoder_params.manifold            = [pars.spiking_inputs{1}(1:end-7) '_pca'];
pars.decoder_params.mani_dims           = pars.mani_dims;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOAD THE DATA
%


% -------------------------------------------------------------------------
% Datafiles for each monkey

files_chewie        = { ...
                        'Chewie_CO_VR_2016-09-12.mat', ...
                        'Chewie_CO_VR_2016-09-14.mat', ...
                        'Chewie_CO_FF_2016-09-15.mat', ...
                        'Chewie_CO_FF_2016-09-19.mat', ...
                        'Chewie_CO_FF_2016-09-21.mat', ...
                        'Chewie_CO_FF_2016-10-05.mat', ...
                        'Chewie_CO_VR_2016-10-06.mat', ...
                        'Chewie_CO_FF_2016-10-07.mat', ...
                        'Chewie_CO_FF_2016-10-11.mat', ...
                        'Chewie_CO_FF_2016-10-13.mat' ...
                        };

files_mihili        = { ...
                        'Mihili_CO_FF_2014-02-03.mat', ...
                        'Mihili_CO_FF_2014-02-17.mat', ...
                        'Mihili_CO_FF_2014-02-18.mat', ...
                        'Mihili_CO_VR_2014-03-03.mat', ...
                        'Mihili_CO_VR_2014-03-04.mat', ...
                        'Mihili_CO_VR_2014-03-06.mat', ...
                        'Mihili_CO_FF_2014-03-07.mat' ...
                        };
            
files_han           = { 'all_TDs_Han.mat' }; 


% -------------------------------------------------------------------------
% Go to the path where the data are, and update the file list if necessary
% (i.e., if using "unsorted" data)

here                = pwd;

% If we want to "unsort", choose the appropriate file names
switch pars.monkey
    % Chewie
    case 'chewie'
        if ~pars.unsorted_yn
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie/Unsorted')
            for f = 1:length(files_chewie), 
                files_chewie{f} = [files_chewie{f}(1:end-4) '_unsorted.mat']; 
            end
        end
        files = files_chewie;
    % Mihili
    case 'mihili'
        if ~pars.unsorted_yn 
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili/Unsorted')
            for f = 1:length(files_mihili), 
                files_mihili{f} = [files_mihili{f}(1:end-4) '_unsorted.mat']; 
            end
        end
        files = files_mihili;
    % Han
    case 'han'
        if ~pars.unsorted_yn
            warning('Data for Han are not sorted'); 
            warning('Loading threshold crossings'); 
            pause(3);
            pars.unsorted_yn = true;
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
            files_han = 'all_TDs_Han.mat';
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
            files_han = 'all_TDs_Han.mat';
        end
        files = files_han;
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
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)} );  
    
    % Keep PCA pars, get rid of the rest of the extra outputs
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end};
    end
    clear pars_td s;

% Unsort the units                        
else      
    
    % % This loads and only does a little bit of preprocessing, with no dim
    % reduction: 1) Loads Baseline trials only (no force field or visuomotor
    % adaptation); 2) Loads only successful trials; 3) Merges all units in
    % the same channel for Mihili and Chewie, but not for Han
    if sum(strcmp(pars.monkey,{'mihili','chewie'}))
        
        master_td 	= loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@mergeUnits});

    % FOR HAN
    else
%         % MERGE UNITS BREAKS IN ONE DATA SET
%         master_td       = loadTDfiles(  files, ...
%                             {@getTDidx,'result','R'}, ...
%                             {@mergeUnits});
        master_td   = loadTDfiles(  files, ...
                            {@getTDidx,'result','R'});
    end

    % Option: only keep electrodes that have units in all the sessions
    if pars.only_common_elecs
        warning('Using only units that are common across sessions'); pause(5);
        master_td   = getCommonUnits(master_td);
    end
end


% go back to where you were path-wise
cd(here);
clear files* here;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET METADATA
%


% GET THE SESSIONS 
meta.sessions       = unique({master_td.date});
n_sessions          = length(meta.sessions);


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


% GET THE TARGETS

% --For Han, we need a couple of tricks
if strcmp(pars.monkey,'han')
    % there are trials with targets = NaN -> exclude them
    master_td       = master_td(~isnan([master_td.target_direction]));
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
end


meta.targets        = unique([master_td.target_direction]);






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (SOME MORE) PRE-PROCESSING
%


% -------------------------------------------------------------------------
% REMOVE CHANNELS WITH LOW FR OR THAT ARE SHUNTED
% -- always done for single neurons, optional for threshold crossings
if pars.unsorted_yn
    pars.remove_bad_neurons
    for s = 1:n_sessions
        this_s      = getTDidx(master_td,{'date',meta.sessions{s}});
        master_td(this_s) = removeBadNeurons(master_td(this_s),pars.bad_neuron_params);
    end
    clear this_s;
end


% -------------------------------------------------------------------------
% FOR THE UNSORTED DATA: 
% 1) SSQUARE-ROOT TRANSFORM AND SMOOTH
% 2) DO PCA 
% (For the spikes we do this during the loading)
if pars.unsorted_yn
    
    % Downsample
    master_td       = binTD(master_td,pars.n_bins_downs);
    
    % Square-root transform
    master_td       = sqrtTransform(master_td,pars.spiking_inputs);
    % Smooth
    master_td       = smoothSignals(master_td,struct(...
                                'signals',pars.spiking_inputs,...
                                'kernelSD',pars.kernel_SD));
    
    
    for s = 1:n_sessions
        this_s = getTDidx( master_td,{'date',meta.sessions{s}} );
        [master_td_new(this_s), pca_info(s)] = getPCA( master_td(this_s), ...
                                    struct(...
                                    'signals',pars.spiking_inputs ));
                                
    end
    master_td       = master_td_new;
%     master_td = binTD(master_td,pars.n_bins_downs);
    clear master_td_new;
end


% -------------------------------------------------------------------------
% Trim the trials

% get the minimum trial duration, if needed for idx_end
if isempty(pars.idx_end{1})
    min_mov_duration = min( arrayfun( @(x) x.idx_trial_end-x.idx_movement_on, master_td) );
    pars.idx_end = {'idx_movement_on',min_mov_duration};
end    

master_td           = trimTD( master_td, pars.idx_start, pars.idx_end );



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EQUALIZE THE NUMBER OF TRIALS ACROSS TARGETS AND SESSIONS
%

master_td           = equalNbrTrialsSessions(master_td);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ESTIMATE DIMENSIONALITY
%

% NEEDTOFIX: SINCE THE BINS ARE PRETTY WIDE AND WE HAVE ONLY 4 TRIALS AND A
% LOT OF NEURONS, THE PCA FUNCTION DOESN'T WORK 
for s = 1:n_sessions
    
    this_s          = getTDidx(master_td,{'date',meta.sessions{s}});
    dims(s)         = estimateDimensionality(master_td(this_s),struct(...
                            'signals',pars.spiking_inputs,...
                            'alpha',0.95,...
                            'num_iter',1000));
                        
    figure, histogram(dims,0:10,'facecolor','k')
    set(gca,'TickDir','out','FontSize',14), box off
    xlabel('Dimensionality'),ylabel('Counts'),
    title([pars.monkey ' - ' pars.spiking_inputs{1}(1:end-7)])
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ALIGN LATENT ACTIVITY OVER SESSIONS
%

aligned_latent_results = align_latent_activity( master_td, pars.align_latent_params );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DECODE KINEMATICS TRAINING BUILDING A DECODER ON ONE DAY AND TESTING IT
% ON ANOTHER
%

decoding_results = decode_from_aligned( master_td, pars.decoder_params );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOTTING
%

SOT_Fig_3_aligned_latent_activity;
