%
% Create and save TD files with threshold crossings rather than sorted
% neurons
%
% ToDo:
%   - Add some fields that Matt was using: perturbation, perturbation_info
%   - See why we can't do removeBadNeurons


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions    

monkey = 'mihili'; % 'mihili'; 'chewie' 


% sessions we want to load
sessions_mihili = { ...
                '2014-02-03', ...
                '2014-02-17', ...
                '2014-02-18', ...
                '2014-03-03', ...
                '2014-03-04', ...
                '2014-03-06', ...
                '2014-03-07' ...
                };
            
sessions_chewie = { ...
                '2016-09-12', ...
                '2016-09-14', ...
                '2016-09-15', ...
                '2016-09-19', ...
                '2016-09-21', ...
                '2016-10-05', ...
                '2016-10-06', ...
                '2016-10-07', ...
                '2016-10-11', ...
                '2016-10-13' ...
                };


% override monkey name and go to path
this_path = pwd;

switch monkey
    case 'chewie'
        data_folder = '/Volumes/I Me Mine/Data/Chewie_8I2';
        sessions = sessions_chewie;
        dest_folder = '/Users/juangallego/Documents/NeuroPlast/Data/Chewie/Unsorted';
    case 'mihili'
        data_folder = '/Volumes/I Me Mine/Data/Mihili_12A3';
        sessions = sessions_mihili;
        dest_folder = '/Users/juangallego/Documents/NeuroPlast/Data/Mihili/Unsorted';
end

cd(data_folder);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now load each file in each folder (session)


% Define parameters
trial_results = {'R','F','I'};
event_list = { ...
    'startTime','trial_start'; ...
    'tgtOnTime','target_on'; ...
    'goCueTime','go_cue'; ...
    'endTime','trial_end'};
td_params = struct( ...
    'event_list',{event_list}, ...
    'trial_results',{trial_results}, ... % which to include
    'exclude_units',[255], ... % sort codes to exclude
    'all_points',true ... % if you want continuous data
);
meta.arrays = {'M1','PMd'};
td_params.meta = meta;


% do
for s = 1:length(sessions)
   
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go to the folder where session 's' is, load each file, and get the TD
    disp(['Loading session ' sessions{s}]);
    cd(sessions{s});
    datafiles = dir(pwd);
    
    for d = 1:length(datafiles)
        
        if datafiles(d).name(1) ~= '.'
            
            load(datafiles(d).name);
            td = parseFileByTrial(cds,td_params);
            
            % add the condition: BL, AD, WO
            for t = 1:length(td)
                td(t).epoch = datafiles(d).name(end-18:end-17);
            end
            
            % STITCH THEM TOGETHER when they are saved
            if exist('trial_data','var')
                trial_data(length(trial_data)+1:length(trial_data)+length(td)) = td;
            else
                trial_data = td;
            end
        end
    end
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADD SOME INFO ABOUT THE TRIAL
    
    % Remove bad trials and reorder fields
    trial_data = removeBadTrials(trial_data);

    % Add movement onset and peak speed
    trial_data = getMoveOnsetAndPeak(trial_data);
    
    trial_data = reorderTDfields(trial_data);
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store data
    file_save = [trial_data(1).monkey '_' trial_data(1).task '_' datafiles(end).name(end-21:end-20) '_',...
        trial_data(1).date(end-3:end) '-' trial_data(1).date(end-9:end-8) '-' trial_data(1).date(end-6:end-5) '_unsorted'];
    
    disp(['Saving ' file_save]);
    save([dest_folder filesep file_save],'trial_data');
    cd(data_folder);
    clear trial_data
end


cd(this_path)