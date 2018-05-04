

root = '/Users/juangallego/Documents/NeuroPlast/Data/Han/';

fd = dir(root);

monkey = 'Han';


% Find TD files in director
idx_td_file = arrayfun(@(x) length(x.name)>3, fd);
fd(~idx_td_file) = [];

non_td_file = ~arrayfun(@(x) strncmp(x.name,monkey,length(monkey)), fd);
fd(non_td_file) = [];


% Load them
this_dir = pwd;
cd(root);

master_td = loadTDfiles({fd.name});


% Get rid of the fields that are not always there (emgs)
master_td = rmfield(master_td,'emg');
master_td = rmfield(master_td,'emg_names');


% Save
trial_data = master_td;
save('all_TDs_Han','trial_data','-v7.3');


cd(this_dir);