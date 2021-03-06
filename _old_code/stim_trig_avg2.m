% 
% Function to do Stimulus Triggered Averaging of an intracortical electrode
% using the Grapevine. 
%
%       function varargout   = stim_trig_avg( varargin )
%
%
% Syntax:
%       EMG                                     = STIM_TRIG_AVG( VARAGIN )
%       [EMG, STA_PARAMS]                       = STIM_TRIG_AVG( VARAGIN )
%       [EMG, FORCE, STA_PARAMS]                = STIM_TRIG_AVG( VARAGIN ),
%               if sta_params.record_force_yn = true
%       [EMG, STA_PARAMS, STA_METRICS]          = STIM_TRIG_AVG( VARAGIN ),
%               if sta_params.record_force_yn = false
%       [EMG, FORCE, STA_PARAMS, STA_METRICS]   = STIM_TRIG_AVG( VARAGIN )
%
%
% Input parameters: 
%       'sta_params'        : stimulation settings. If not passed, read
%                               from stim_trig_avg_default 
% Outputs: 
%       'emg'               : EMG evoked by each simulus, and some
%                               related information
%       'force'             : Force evoked by each stimulus, similar to
%                               the emg field. This is the second para
%       'sta_params'        : stimulation parameters
%       'sta_metrics'       : StTA metrics computed with the
%       'calculate_sta_metrics' function  
%
%
%
%                           Last modified by Juan Gallego 6/17/2015




% %%%%%%%%%%%%
%   ToDo: 
%       - When the first threshold crossing is missed, read the rest of the
%       data
%   Some known issues:
%       - The resistor used in the analog input is hard-coded (100 Ohm)





function varargout = stim_trig_avg2( varargin )


close all;


% read parameters

if nargin > 1   
    error('ERROR: The function only takes one argument of type StTA_params');
elseif nargin == 1
    sta_params                  = varargin{1};
elseif nargin == 0
    sta_params                  = stim_trig_avg_default2();
end

if nargout > 4
    disp('ERROR: The function only returns up to three variables, of type emg, force and sta_params');
end



%--------------------------------------------------------------------------
%% connect with Central 

% connect to central; if connection fails, return error message and quit
if ~cbmex('open', 1)
    
    echoudp('off');
%    close(handles.keep_running);
    error('ERROR: Connection to Central Failed');
end



% If we want to save the data ...

% Note structure 'hw' will have all the cerebus and grapevine stuff
if sta_params.save_data_yn
   
    % create file name
    hw.data_dir                 = [sta_params.data_dir filesep 'STA_data_' datestr(now,'yyyy_mm_dd')];
    if ~isdir(hw.data_dir)
        mkdir(hw.data_dir);
    end
    hw.start_t                  = datestr(now,'yyyymmdd_HHMMSS');
    hw.cb.full_file_name        = fullfile( hw.data_dir, [sta_params.monkey '_' sta_params.bank '_' num2str(sta_params.stim_elec) '_' hw.start_t '_' sta_params.task '_STA' ]);

    % start 'file storage' app, or stop ongoing recordings
    cbmex('fileconfig', fullfile( hw.data_dir, hw.cb.full_file_name ), '', 0 );  
    drawnow;                        % wait till the app opens
    pause(1);
    drawnow;                        % wait some more to be sure. If app was closed, it did not always start recording otherwise

    % start cerebus file recording
    cbmex('fileconfig', hw.cb.full_file_name, '', 1);
end


% configure acquisition
cbmex('trialconfig', 1);            % start data collection
drawnow;

pause(1);                           % ToDo: see if it's necessary



% figure out how many EMG channels there are, and preallocate matrices
% accordingly. Check that there is a 'sync out' signal from the stimulator 

[ts_cell_array, ~, analog_data] = cbmex('trialdata',1);

% look for the 'sync out' signal ('Stim_trig')
hw.cb.stim_trig_ch_nbr  =  find(strncmp(ts_cell_array(:,1),'Stim',4));
if isempty(hw.cb.stim_trig_ch_nbr)
    error('ERROR: Sync signal not found in Cerebus. The channel has to be named Stim_trig');
else
    disp('Sync signal found');
    
    % define resistor to record sync pulse
    hw.cb.sync_out_resistor = 100;
end



% EMG data will be stored in the 'emg' data structure
% 'emg.evoked_emg' has dimensions EMG signal -by- EMG channel- by- stimulus
% nbr 

analog_data(:,1)            = ts_cell_array([analog_data{:,1}]',1); % ToDo: replace channel numbers with names
emg.labels                  = analog_data( strncmp(analog_data(:,1), 'EMG', 3), 1 );
emg.nbr_emgs                = numel(emg.labels); disp(['Nbr EMGs: ' num2str(emg.nbr_emgs)]), disp(' ');
emg.fs                      = cell2mat(analog_data(find(strncmp(analog_data(:,1), 'EMG', 3),1),2));
emg.length_evoked_emg       = ( sta_params.t_before + sta_params.t_after ) * emg.fs/1000 + 1;
emg.evoked_emg              = zeros( emg.length_evoked_emg, emg.nbr_emgs, sta_params.nbr_stims_ch ); 



% If chosen to record force

if sta_params.record_force_yn
   
    force.labels            = analog_data( strncmp(analog_data(:,1), 'Force', 5), 1 );
    force.nbr_forces        = numel(force.labels); disp(['Nbr Force Sensors: ' num2str(force.nbr_forces)]), disp(' ');
    force.fs                = cell2mat(analog_data(find(strncmp(analog_data(:,1), 'Force', 5),1),2));
    force.length_evoked_force   = ( sta_params.t_before + sta_params.t_after ) * force.fs/1000 + 1;
    force.evoked_force      = zeros( force.length_evoked_force, force.nbr_forces, sta_params.nbr_stims_ch );
end


clear analog_data ts_cell_array;
cbmex('trialconfig', 0);        % stop data collection until the stim starts



%--------------------------------------------------------------------------
%% connect with Grapevine

% initialize xippmex
hw.gv.connection            = xippmex;

if hw.gv.connection ~= 1
    cbmex('close');
    error('ERROR: Xippmex did not initialize');
end


% check if the sync out channel has been mistakenly chosen for stimulation
if ~isempty(find(sta_params.stim_elec == sta_params.sync_out_elec,1))
    cbmex('close');
    error('ERROR: sync out channel chosen for ICMS!');
end


% find all Micro+Stim channels (stimulation electrodes). Quit if no
% stimulator is found 
hw.gv.stim_ch               = xippmex('elec','stim');

if isempty(hw.gv.stim_ch)
    cbmex('close');
    error('ERROR: no stimulator found!');
end


% quit if the specified channel (in 'sta_params.stim_elec') does not exist,
% or if the sync_out channel does not exist  

if isempty(find(hw.gv.stim_ch==sta_params.stim_elec,1))
    cbmex('close');
    error('ERROR: stimulation channel not found!');
elseif isempty(find(hw.gv.stim_ch==sta_params.sync_out_elec,1))
    cbmex('close');
    error('ERROR: sync out channel not found!');
end


% SAFETY! check that the stimulation amplitude is not too large ( > 90 uA
% or > 1 ms) 
if sta_params.stim_ampl > 0.090
    cbmex('close');
    error('ERROR: stimulation amplitude is too large (> 90uA) !');    
elseif sta_params.stim_pw > 1
    cbmex('close');
    error('ERROR: stimulation pulse width is too large (> 1ms) !');    
end
   


%--------------------------------------------------------------------------
%% some preliminary stuff


% this creates 'epochs' of continuous ICMS and cerebus recordings. To avoid
% the need of reading a huge chunk of data from Central at once
hw.cb.epoch_duration        = 10;   % epoch duration (in s)
hw.cb.nbr_epochs            = ceil(sta_params.nbr_stims_ch/sta_params.stim_freq/hw.cb.epoch_duration);
hw.cb.nbr_stims_this_epoch  = sta_params.stim_freq*hw.cb.epoch_duration;
hw.cb.ind_ev_emg            = 0;    % ptr to know where to store the evoked EMG


drawnow;


%--------------------------------------------------------------------------
%% stimulate to get STAs



%------------------------------------------------------------------
% Define the stimulation string and start data collection
% Note that TD adds a delay that is the time before the stimulation for
% the StTA (defined in sta_params.t_before) + 10 ms, to avoid
% synchronization issues

stim_string             = [ 'Elect = ' num2str(sta_params.stim_elec) ',' num2str(sta_params.sync_out_elec) ',;' ...
                            'TL = ' num2str(1000/sta_params.stim_freq) ',' num2str(1000/sta_params.stim_freq) ',; ' ...
                            'Freq = ' num2str(sta_params.stim_freq) ',' num2str(sta_params.stim_freq) ',; ' ...
                            'Dur = ' num2str(sta_params.stim_pw) ',' num2str(sta_params.stim_pw) ',; ' ...
                            'Amp = ' num2str(sta_params.stim_ampl/sta_params.stimulator_resolut) ',' num2str(ceil(3/hw.cb.sync_out_resistor/sta_params.stimulator_resolut*1000)) ',; ' ...
                            'TD = ' num2str((sta_params.t_before)/1000) ',' num2str((sta_params.t_before)/1000) ',; ' ...
                            'FS = 0,0,; ' ...
                            'PL = 1,1,;'];


for i = 1:hw.cb.nbr_epochs

    % start data collection
    cbmex('trialconfig', 1);
    drawnow;
    drawnow;
    drawnow;
    
    %------------------------------------------------------------------
    % Stimulate the channel as many times as specified

    if i == hw.cb.nbr_epochs 
        hw.cb.nbr_stims_this_epoch  = rem(sta_params.nbr_stims_ch,sta_params.stim_freq*30);
    end

    for ii = 1:hw.cb.nbr_stims_this_epoch

        t_start             = tic;
        drawnow;

        % send stimulation command
        xippmex('stim',stim_string);
        drawnow;
        drawnow;
        drawnow;


        % wait for the inter-stimulus interfal (defined by stim freq)
        t_stop              = toc(t_start);
        while t_stop < (1/sta_params.stim_freq)
            t_stop          = toc(t_start);
        end

        % if this is the last stimulus in this epoch, record for some
        % extra time. To avoid loosing the last epoch ?it may happen
        % some times  
        if ii == hw.cb.nbr_stims_this_epoch
            pause(0.1);
        end
    end


    %------------------------------------------------------------------
    % read EMG data and sync pulses

    % read the data from central (flush the data cache)
    [ts_cell_array, ~, analog_data] = cbmex('trialdata',1);
    cbmex('trialconfig', 0);
    drawnow;

    % display if some of the stim pulses were lost
    if numel(cell2mat(ts_cell_array(hw.cb.stim_trig_ch_nbr,2))) ~= hw.cb.nbr_stims_this_epoch
        disp(' ');
        disp(['Warning: ' num2str(hw.cb.nbr_stims_this_epoch - numel(cell2mat(ts_cell_array(hw.cb.stim_trig_ch_nbr,2)))) ' sync pulses lost in this stim epoch!']);
    else
        disp(' ');
        disp('all sync pulses detected for this stim epoch');
    end

    
    %------------------------------------------------------------------
    % retrieve EMG (and Force, if specified) data, and stimulation time
    % stamps 
    analog_data(:,1)        = ts_cell_array([analog_data{:,1}]',1);
    aux                     = analog_data( strncmp(analog_data(:,1), 'EMG', 3), 3 );

    for ii = 1:emg.nbr_emgs
        emg.data(:,ii)   = double(aux{ii,1}); 
    end
    clear aux
    
    if sta_params.record_force_yn
        aux2                = analog_data( strncmp(analog_data(:,1), 'Force', 5), 3 ); % ToDo: double check this line
        for ii = 1:force.nbr_forces
            force.data(:,ii)   = double(aux2{ii,1});
        end
    end

    ts_sync_pulses          = double(cell2mat(ts_cell_array(hw.cb.stim_trig_ch_nbr,2)));

    
    
    % RESOLVE ISSUES RELATED TO THE SYNCHRONIZATION BETWEEN TIME STAMPS AND
    % ANALOG SIGNALS WHEN READING FROM CENTRAL. THIS HAS BEEN FIXED IN
    % CBMEX v6.3, ALTHOUGH IT SOMETIMES MISSES THE FIRST THRESHOLD CROSSING  
    
    ts_sync_pulses_its_freq             = ts_sync_pulses/30000*emg.fs;
    %    ts_first_sync_pulse_emg_freq        = ts_sync_pulses_its_freq(find(ts_sync_pulses_its_freq<5000,1));
    analog_sync_signal                  = double(analog_data{10,3});
    
    ts_first_sync_pulse_analog_signal   = find( (analog_sync_signal-mean(analog_sync_signal)) <-1.5*std(analog_sync_signal), 1);


    if abs( ts_first_sync_pulse_analog_signal - ts_sync_pulses_its_freq(1) ) > emg.fs/1000
        if ts_first_sync_pulse_analog_signal < ts_sync_pulses_its_freq(1)
            disp('Central has skipped the first threshold crossing of the sync signal!!!');
            disp(['the delay btw analog/thr crossing is: ' num2str( (ts_first_sync_pulse_analog_signal - ts_sync_pulses_its_freq(1))/10 )])
        else
            disp('the delay between the time stamps and the analog signal is > 1 ms!!!');
            disp(['it is: ' num2str( (ts_first_sync_pulse_analog_signal - ts_sync_pulses_its_freq(1))/10 )])
        end
    else
        disp('the delay between the time stamps and the analog signal is < 1 ms');
    end

    
%     figure,plot(analog_sync_signal), hold on, xlim([0 10000]), xlabel(['sample numer at EMG fs = ' num2str(emg.fs) ' (Hz)']), 
%     stem(ts_sync_pulses_its_freq,ones(length(ts_sync_pulses),1)*-5000,'marker','none','color','r'), legend('analog signal','time stamps')
%    ToDo: DELETE UNTIL HERE


    % when the time stamps and the analog data are not synchronized, they
    % won't be stored; 
    if abs( ts_first_sync_pulse_analog_signal - ts_sync_pulses_its_freq(1) ) <  emg.fs/1000


        % remove the first and/or the last sync pulse if they fall outside the data
        if floor(ts_sync_pulses(1)/30000) < sta_params.t_before/1000
            ts_sync_pulses(1)   = [];
        elseif (length(emg.data) - sta_params.t_after) < floor(ts_sync_pulses(end)/30000*emg.fs)
            ts_sync_pulses(end) = [];
        end


        %------------------------------------------------------------------
        % store the evoked EMG (interval around the stimulus defined by
        % t_before and t_after in params 

        for ii = 1:min(length(ts_sync_pulses),length(emg.evoked_emg))
            trig_time_in_emg_sample_nbr     = floor(double(ts_sync_pulses(ii))/30000*emg.fs - sta_params.t_before/1000*emg.fs);

            % check if we haven't recorded EMG for long enough, during some of
            % the stimuli. If not, store the data
            if (trig_time_in_emg_sample_nbr + (sta_params.t_after + sta_params.t_before)*emg.fs/1000 ) < length(emg.data)
                emg.evoked_emg(:,:,ii+hw.cb.ind_ev_emg)    = emg.data( trig_time_in_emg_sample_nbr : ...
                    (trig_time_in_emg_sample_nbr + emg.length_evoked_emg - 1), : );
            else
                disp('one sync pulse in the EMG is far too late!');
                drawnow;
            end
        end

        
        %------------------------------------------------------------------
        % store the evoked Force (interval around the stimulus defined by
        % t_before and t_after in params
        
        if sta_params.record_force_yn
            
            for ii = 1:min(length(ts_sync_pulses),length(emg.evoked_emg))
                
                trig_time_in_force_sample_nbr   = floor(double(ts_sync_pulses(ii))/30000*force.fs - sta_params.t_before/1000*force.fs);
                
                % check if we haven't recorded Force for long enough, during some
                % of the stimuli. If not, store the data
                if (trig_time_in_force_sample_nbr + (sta_params.t_after + sta_params.t_before)*force.fs/1000 ) < length(force.data)
                    force.evoked_force(:,:,ii+hw.cb.ind_ev_emg)   = force.data( trig_time_in_force_sample_nbr : ...
                        (trig_time_in_force_sample_nbr + force.length_evoked_force - 1), : );
                else
                    disp('one sync pulse in the Force is far too late!');
                    drawnow;
                end
            end
        end
    end

    
    
    % update ptr to index
    hw.cb.ind_ev_emg        = length(ts_sync_pulses) + hw.cb.ind_ev_emg;


    % delete some variables
    clear analog_data ts_cell_array; 
    emg                     = rmfield(emg,'data');
    if sta_params.record_force_yn
        force                   = rmfield(force,'data');
    end
end




%--------------------------------------------------------------------------
% Save data and stop cerebus recordings


disp(['Finished stimulating electrode ' num2str(sta_params.stim_elec)]);
disp(' ');


% Save the data, if specified in sta_params
if sta_params.save_data_yn
    
    % stop cerebus recordings
    cbmex('fileconfig', hw.cb.full_file_name, '', 0);
    cbmex('close');
    drawnow;
    drawnow;
    disp('Communication with Central closed');
    
%    xippmex('close');

    % save matlab data. Note: the time in the faile name will be the same as in the cb file
    hw.matlab_full_file_name    = fullfile( hw.data_dir, [sta_params.monkey '_' sta_params.bank '_' num2str(sta_params.stim_elec) '_' hw.start_t '_' sta_params.task '_STA' ]);
    
    disp(' ');
    
    if sta_params.record_force_yn == false
        save(hw.matlab_full_file_name,'emg','sta_params');
        disp(['EMG data and Stim Params saved in ' hw.matlab_full_file_name]);
    else
        save(hw.matlab_full_file_name,'emg','force','sta_params');
        disp(['EMG and Force data and Stim Params saved in ' hw.matlab_full_file_name]);
    end    
end

cbmex('close')


% Calculate the STA metrics and plot, if specified in sta_params
if sta_params.plot_yn
   
    if sta_params.record_force_yn == false
        sta_metrics             = calculate_sta_metrics2( emg, sta_params );
    else
        sta_metrics             = calculate_sta_metrics2( emg, force, sta_params );
    end
end



%-------------------------------------------------------------------------- 
% Return variables
if nargout == 1
    varargout{1}            = emg;
elseif nargout == 2
    varargout{1}            = emg;
    varargout{2}            = sta_params;
elseif nargout == 3
    if sta_params.record_force_yn     
        varargout{1}        = emg;
        varargout{2}        = force;
        varargout{3}        = sta_params;
    else
        varargout{1}        = emg;
        varargout{2}        = sta_params;
        varargout{3}        = sta_metrics;        
    end
elseif nargout == 4
    varargout{1}            = emg;
    varargout{2}            = sta_params;
    varargout{3}            = emg;
    varargout{4}            = sta_metrics;
end


