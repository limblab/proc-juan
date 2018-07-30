%
% Rename fields in td structs so it follows Matt's convention
%


clear, close all; 


% will load all the files in the path if you don't pass the files
path = '/Users/juangallego/Documents/NeuroPlast/Data/Han'; %'/Users/juangallego/Documents/NeuroPlast/Data/Chips';

save_yn = true;


% Raeed's names
raeed_names = { 'idx_startTime', ...
                'idx_goCueTime', ...
                'idx_tgtOnTime', ...
                'idx_endTime', ...
                'idx_bumpTime', ...
                };
            
% Standard names
std_names   = { 'idx_trial_start', ...
                'idx_go_cue', ...
                'idx_target_on', ...
                'idx_trial_end', ...
                'idx_bump' ...
                };
                
            
            
disp(['Renaming the IDX field of all the files in ' path]);

here = pwd; cd(path);
fn = dir;

% delete files with very short names --can't be m-files
fn(arrayfun(@(x) length(x.name)<4, fn)) = [];
% find m_files
mfn = fn( arrayfun(@(x) strcmp(x.name(end-3:end),'.mat'), fn) );


% Load all the M-files in the folder
for m = 1:length(mfn)

    load(mfn(m).name);
    disp(['Fixing ' mfn(m).name]);

    allfields = fieldnames(trial_data);


    % Loop through all the names that need to be fixed
    for r = 1:length(raeed_names)

        % find this field in the data
        idx_tfield = strncmp(allfields,raeed_names{r},length(raeed_names{r}));

        if sum(idx_tfield) == 1

            % retrieve the standard name for the list
            std_name = std_names{r};

            % create field with standard name and the right data
            for t = 1:length(trial_data)

                trial_data(t).(std_name) = trial_data(t).(allfields{idx_tfield});
            end

            % remove old named field
            trial_data = rmfield(trial_data,allfields{idx_tfield});
            trial_data = reorderTDfields(trial_data);
        end
    end

    % Save?
    if save_yn
        disp(['Saving ' mfn(m).name]);
        save(mfn(m).name,'trial_data','-v7.3');
    end
end
