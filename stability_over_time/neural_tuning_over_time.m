clear;
clc;
close all;

% simple tuning script
pars.bad_neuron_params.min_fr           = 0.1;
% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;

pars.tuning_window = {'idx_go_cue','idx_trial_end'};

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

filedir = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/';

files = files_mihili;

filepaths = cell(size(files));
for iFile = 1:length(files)
    filepaths{iFile} = fullfile(filedir,files{iFile});
end


trial_data = loadTDfiles(filepaths, ...
    {@getTDidx,'epoch','BL','result','R'}, ...
    @stripSpikeSorting, ...
    {@removeBadNeurons,pars.bad_neuron_params}, ...
    {@trimTD,pars.tuning_window{1},pars.tuning_window{2}});

%% do the tuning
file_sgs = cell(length(files),1);
file_tuning = cell(length(files),3);
for iFile = 1:length(files)
    % get the data for this file
    [~,temp,~] = fileparts(files{iFile});
    temp = strsplit(temp,'_');
    the_date = temp{4};
    
    [~,td] = getTDidx(trial_data,'date',datestr(the_date,'mm-dd-yyyy'));
    
    % get average firing rate in window
    fr = cell2mat(cellfun(@(x) sum(x,1)/(0.01*size(x,1)),{td.M1_spikes},'uni',0)');
    % get target directions
    dir = [td.target_direction];
    % fit tuning curves
    [tc,cb,r] = regressTuningCurves(fr,dir',{'bootstrap',100,0.95});
    % store it for later
    file_tuning{iFile,1} = tc;
    file_tuning{iFile,2} = cb;
    file_tuning{iFile,3} = r;
    file_sgs{iFile} = td(1).M1_unit_guide;
    
    
end


%% plot the stuff

figure;
hold all;

for iFile = 1:size(file_tuning,1)
    [elecs,idx1,idx2] = intersect(file_sgs{1}(:,1),file_sgs{iFile}(:,1));

    d = angleDiff(file_tuning{1,1}(idx1,3),file_tuning{iFile,1}(idx2,3),true,false);
    
    m = 180/pi*circular_mean(d);
    s = 180/pi*circular_std(d)/sqrt(length(d));
    plot(iFile,m,'ko','LineWidth',2);
    plot([iFile, iFile],[m-s,m+s],'k-','LineWidth',2)
end
set(gca,'Box','off','TickDir','out','FontSize',14);
title('Mihili');
xlabel('Session Number')
ylabel('| Change in Preferred Direction Compared to Day 1 |');



