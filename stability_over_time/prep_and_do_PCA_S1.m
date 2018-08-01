%
% Finish processing S1 monkeys -- prepares the data and does PCA for the S1
% monkeys, because we can't do it in the loader due to some issues with the
% TD files
%


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