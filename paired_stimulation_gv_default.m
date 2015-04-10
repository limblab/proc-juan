%
% Parameters for paired intracortical-muscle stimulation using
% paired_stimulation() 
%
%       function p_s_params = paired_stimulation_default()
%
%       p_s_params              : structure with fields
%           ISI                 : inter-stimulus interval. + = cortex first (ms)
%           t_btw_pairs         : time btw pairs of stimulus (wrt the first stimulus in the pair; in ms)
%           cortical_elec       : cortical electrode for ICMS
%           stim_ampl_cx        : stim amplitude for Cx (mA)
%           stim_freq_cx        : stim frequency for Cx (Hz)
%           stim_pw_cx          : stim pulse width for Cx [and sync out] (ms)
%           train_duration_cx   : duration of the train delivered to Cx (ms)
%           muscle_elec         : muscle electrodes
%           stim_ampl_muscle    : stim amplitude for the muscles (mA)
%           stim_freq_muscle    : stim frequency for the muscles (Hz)
%           stim_pw_muscle      : stim pulse width for the muscle (ms)
%           train_duration_muscle   : duration of the train delivered to the muscles(ms)
%           sync_out_elec       : electrode that generates the sync signal for STA the force data
%           nbr_stimuli         : nbr of paired stimuli
%           pre_stim_win        : time (ms), before the delivery of a stimulus, during which force will be recorded and plotted
%           post_stim_win       : time (ms), after the delivery of a stimulus, during which force will be recorded and plotted
%           sync_out_elec       : ch nbr of the sync out electrode
%           stimulator_resolut  : resolution of the Cx stimulator (mA)
%           record_force_yn     : (=1) record force; (=0) don't
%           save_data_yn        : save Matlab and cerebus data
%           data_dir            : directory where the data will be saved
%           monkey              : monkey name, to generate the filename
%           bank                : bank of the array



function p_s_params = paired_stimulation_gv_default()


p_s_params      = struct(...
    'ISI',                      15, ...
    't_btw_pairs',              100, ...
    'cortical_elec',            15, ...
    'stim_ampl_cx',             0.054, ...
    'stim_freq_cx',             20, ...
    'stim_pw_cx',               0.2, ...
    'train_duration_cx',        50, ...
    'muscle_elec',              [1 2], ...
    'stim_ampl_muscle',         [2 2], ...
    'stim_freq_muscle',         [20 20], ...
    'stim_pw_muscle',           [0.2 0.2], ...
    'stim_polar_muscle',        [0 0], ...
    'train_durration_muscle',   [50 50], ...
    'nbr_stimuli',              1200, ...
    'pre_stim_win',             300, ...
    'post_stim_win',            500, ...
    'sync_out_elec',            32, ...
    'stimulator_resolut',       0.018, ...
    'record_force_yn',          0, ...
    'save_data_yn',             0, ...
    'data_dir',                 'c:\Users\limblab\Desktop\temp code', ...
    'monkey',                   'Jango', ...
    'bank',                     'A');