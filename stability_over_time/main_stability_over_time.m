
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

pars.monkey         = 'Chewie'; % 'chewie2'; 'chewie'; 'mihili'; 'han'; 'chips'; 'jaco'
pars.spiking_inputs = {'PMd_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}; {'S1_spikes'}

% Sesssions to discard if any
pars.sessions_discard = []; %6:14; %[12 13 14];

% -------------------------------------------------------------------------
% Load rest of the default parameters
params_stability_over_time;


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
        files = files_chewie;
    % Chewie M1 only implant
    case 'chewie2'
        files = files_chewie2;
    case 'jaco'
        files = files_jaco;
    % Mihili M1-PMd
    case 'mihili'
        files = files_mihili;
    % MrT PMd
    case 'mrt'
        files = files_mrt;
    % Han S1
    case 'han'
        if ~pars.unsorted_yn
            warning('Data for Han are not sorted; Loading threshold crossings'); pars.unsorted_yn = true;
        end
        files = files_han;
    % Chips S1
    case 'chips'
        if ~pars.unsorted_yn
            warning('Data for Chips are not sorted; Loading threshold crossings'); pars.unsorted_yn = true;
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
    
    if sum(strcmpi(pars.monkey,{'mihili','chewie','chewie2','mrt','jaco'}))

        % FOR PMd
        if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
        
            [master_td, pars_td] = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL','result','R'}, ...
                            {@stripSpikeSorting},...
                            {@removeBadTrials,struct('nan_idx_names','idx_movement_on')}, ...
                            {@binTD,pars.n_bins_downs}, ...
                            {@removeBadNeurons,pars.bad_neuron_params},...
                            {@sqrtTransform,pars.spiking_inputs}, ...
                            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
                            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
                            {@getPCA,struct('signals',pars.spiking_inputs)});
        master_td = rmfield(master_td,'force');
        % FOR M1                        
        elseif strcmpi(pars.spiking_inputs{1},'M1_spikes')
            
            % For Chewie and Mihili things work without removing any trials
            if strcmpi(pars.monkey,'chewie') || strcmpi(pars.monkey,'mihili')
                
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
                        
                 % some sessions are missing force so just throw it out
                 % since we don't n eed it for this analysis
                 master_td = rmfield(master_td,'force');
            % For Chewie2 and Jaco we need to remove some trials 
            elseif strcmpi(pars.monkey,'chewie2') || strcmpi(pars.monkey,'jaco')
         
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
                            {@getPCA,struct('signals',pars.spiking_inputs)}, ...
                            {@trimTD, pars.idx_start, pars.idx_end} );
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
if strcmpi(pars.spiking_inputs{1},'PMd_spikes')
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end};
    end
elseif ~( strcmpi(pars.monkey,'han') || strcmpi(pars.monkey,'chips') )
    for s = 1:size(pars_td.extra_outs,1)
        pca_info(s) = pars_td.extra_outs{s,end-1};
    end
end
clear pars_td s;


% go back to where you were path-wise
cd(here);
clear files* here min_* i;


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
clear i_sort s sorted_dates;


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
    align_results       = align_latent_activity( trimTD(master_td,pars) , pars.align_latent_params );
    within_day_align_results = align_latent_activity_within_day( trimTD(master_td_all_trials,pars), ...
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
    clas_results = classify_across_days(master_td,'align',pars.class_params);

    % now do spikes
    clas_spike_results = classify_across_days(master_td,'spikes',pars.class_params);
    
    % save some results
    fn1 = [pars.monkey '_' pars.spiking_inputs{1}(1:end-7) '_Classification_over_time_' num2str(length(pars.mani_dims)) 'D'];
    
    savefig(f1,fullfile(params.save_dir,pars.spiking_inputs{1}(1:end-7),[fn1 '.fig']));
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

