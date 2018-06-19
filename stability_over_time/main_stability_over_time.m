
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability of latent activity over time --main function

clear all, close all



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%


% -------------------------------------------------------------------------
% What data to use

pars.monkey         = 'chewie'; % 'chewie'; 'mihili'
pars.spiking_inputs = {'M1_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}

% Sesssions to discard if any
pars.sessions_discard = [];

% -------------------------------------------------------------------------
% Data preprocessing

% "Unsort" the sorted neurons?
pars.unsorted_yn    = false; 

% Minimum FR per channel
pars.min_fr         = 0.1;
% Gaussian kernel for smoothing
pars.kernel_SD      = 0.1;

% Neural modes to use
pars.mani_dims      = 1:10; % 1:10; 'estimate';

% Window start & end --if idx_end is empty: the duration of the shortest trial 
pars.idx_start      = {'idx_go_cue',0}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
pars.idx_end        = {'idx_go_cue',18}; % {''}; % {'idx_go_cue',10}

% "Downsampling rate": nbr of bins that will be combined
pars.n_bins_downs   = 5;

% -------------------------------------------------------------------------
% Decoder settings -to predict behavior

% Inputs and Outputs
pars.decoder_out    = 'vel';
pars.decoder_in     = 'aligned_data'; % 'aligned_data'; 'raw_data'

% Bins for decoder history
pars.hist_bins      = 3;

% Folds for multi-fold cross-validation
pars.n_folds        = 6; 



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
switch monkey
    % Chewie
    case 'chewie'
        if ~unsorted_yn
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
        if ~unsorted_yn 
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
        if ~unsorted_yn
            warning('Data for Han are not sorted'); pause(3);
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Han')
        end
end


% -------------------------------------------------------------------------
% Load the data

% Sorted data
if ~unsorted_yn
    
    % This loads and already preprocesses and does dimensionality reduction
    % 1) Loads Baseline trials only (no force field or visuomotor
    % adaptation); 2) Loads only successful trials; 3) Square
    % root-transforms, smooths and does PCA (on the entire trial); 4)
    % Downsamples the data 
    master_td       = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@sqrtTransform,spiking_inputs}, ...
                            {@smoothSignals,struct('signals',spiking_inputs,'calc_fr',true,'kernel_SD',kernel_SD)}, ...
                            {@getPCA,struct('signals',spiking_inputs)}, ...
                            {@binTD,n_bins} ...
                            );    
% Unsort the units                        
else      
    
    % % This loads and only does a little bit of preprocessing, with no dim
    % reduction: 1) Loads Baseline trials only (no force field or visuomotor
    % adaptation); 2) Loads only successful trials; 3) Merges all units in
    % the same channel
    master_td       = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@mergeUnits});
                        
%     master_td   = loadTDfiles(  files, ...
%                             {@getTDidx,'epoch','BL'}, ...
%                             {@getTDidx,'result','R'}, ...
%                             {@mergeUnits},...
%                             {@sqrtTransform,spiking_inputs}, ...
%                             {@smoothSignals,struct('signals',spiking_inputs,'calc_fr',true,'kernel_SD',kernel_SD)}, ...
%                             {@getPCA,struct('signals',spiking_inputs)}, ...
%                             {@binTD,n_bins} ...
%                             );

    % Necessary???
    warning('Using only units that are common across sessions'); pause(5);
    master_td       = getCommonUnits(master_td);
end


% go back to where you were path-wise
cd(here);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (SOME MORE) PRE-PROCESSING
%





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET METADATA
%


% Get the sessions
meta.sessions       = unique({master_td.date});
% Get the targets
meta.targets        = unique({master_td.target_direction});

