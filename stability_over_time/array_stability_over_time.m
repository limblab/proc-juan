clear;
clc;
close all;

% data_root = '/Volumes/MattData/';
% data_save_dir = '/Users/mattperich/Dropbox/Research/Data/StabilityData';
% out_dir = '/Users/mattperich/Dropbox/Research/Papers/2018 - Stability latent activity/Results/Neuron stability/';

data_root = 'C:\Users\Matt\Desktop\LimblabData\';
data_save_dir = 'C:\Users\Matt\Desktop\LimblabData\';
out_dir = 'C:\Users\Matt\Desktop\LimblabData\';

do_data_processing = true; % can only do this if you have raw data
do_tracking        = true; % can only do this if you have completed data processing
save_results       = true; % save a .mat file with the results
save_figs          = true;

pars.monkey = 'Chewie2';
pars.array = 'M1';
files_chewie = { ...
    'Chewie_CO_VR_2016-09-09.mat', ...
    'Chewie_CO_VR_2016-09-12.mat', ...
    'Chewie_CO_VR_2016-09-14.mat', ...
    'Chewie_CO_FF_2016-09-15.mat', ...
    'Chewie_CO_FF_2016-09-19.mat', ...
    'Chewie_CO_FF_2016-09-21.mat', ...
    'Chewie_CO_FF_2016-09-23.mat', ...
    'Chewie_CO_VR_2016-09-29.mat', ...
    'Chewie_CO_FF_2016-10-05.mat', ...
    'Chewie_CO_VR_2016-10-06.mat', ...
    'Chewie_CO_FF_2016-10-07.mat', ...
    'Chewie_CO_FF_2016-10-11.mat', ...
    'Chewie_CO_FF_2016-10-13.mat' ...
    'Chewie_CO_CS_2016-10-14.mat', ...
    'Chewie_CO_CS_2016-10-21.mat', ...
    };

files_chewie2       = { ...
    'Chewie_CO_VR_2013-10-03.mat', ...
    'Chewie_CO_FF_2013-10-22.mat', ...
    'Chewie_CO_FF_2013-10-23.mat', ...
    'Chewie_CO_FF_2013-10-31.mat', ...
    'Chewie_CO_FF_2013-11-01.mat', ...
    'Chewie_CO_FF_2013-12-03.mat', ...
    'Chewie_CO_FF_2013-12-04.mat', ...
    'Chewie_CO_VR_2013-12-19.mat', ...
    'Chewie_CO_VR_2013-12-20.mat', ...
    'Chewie_CO_CS_2015-03-09.mat', ...
    'Chewie_CO_CS_2015-03-11.mat', ...
    'Chewie_CO_CS_2015-03-12.mat', ...
    'Chewie_CO_CS_2015-03-13.mat', ...
    'Chewie_CO_CS_2015-03-19.mat', ...
    'Chewie_CO_FF_2015-06-29.mat', ...
    'Chewie_CO_FF_2015-06-30.mat', ...
    'Chewie_CO_FF_2015-07-01.mat', ...
    'Chewie_CO_FF_2015-07-03.mat', ...
    'Chewie_CO_FF_2015-07-06.mat', ...
    'Chewie_CO_FF_2015-07-07.mat', ...
    'Chewie_CO_FF_2015-07-08.mat', ...
    'Chewie_CO_VR_2015-07-09.mat', ...
    'Chewie_CO_VR_2015-07-10.mat', ...
    'Chewie_CO_VR_2015-07-13.mat', ...
    'Chewie_CO_VR_2015-07-14.mat', ...
    'Chewie_CO_VR_2015-07-15.mat', ...
    'Chewie_CO_VR_2015-07-16.mat', ...
    'Chewie_CO_CS_2015-11-03.mat', ...
    'Chewie_CO_CS_2015-11-04.mat', ...
    'Chewie_CO_CS_2015-11-06.mat', ...
    'Chewie_CO_GR_2015-11-09.mat', ...
    'Chewie_CO_GR_2015-11-10.mat', ...
    'Chewie_CO_GR_2015-11-12.mat', ...
    'Chewie_CO_GR_2015-11-13.mat', ...
    'Chewie_CO_GR_2015-11-16.mat', ...
    'Chewie_CO_GR_2015-11-17.mat', ...
    'Chewie_CO_VR_2015-11-19.mat', ...
    'Chewie_CO_GR_2015-11-20.mat', ...
    'Chewie_CO_VR_2015-12-01.mat', ...
    'Chewie_CO_VR_2015-12-03.mat', ...
    'Chewie_CO_VR_2015-12-04.mat' ...
    };

files_mihili = { ...
    'Mihili_CO_FF_2014-02-03.mat', ...
    'Mihili_CO_FF_2014-02-17.mat', ...
    'Mihili_CO_FF_2014-02-18.mat', ...
    'Mihili_CO_VR_2014-03-03.mat', ...
    'Mihili_CO_VR_2014-03-04.mat', ...
    'Mihili_CO_VR_2014-03-06.mat', ...
    'Mihili_CO_FF_2014-03-07.mat', ...
    'Mihili_CO_CS_2014-09-29.mat', ...
    'Mihili_CO_CS_2014-12-03.mat', ...
    'Mihili_CO_CS_2015-05-11.mat', ...
    'Mihili_CO_FF_2015-06-10.mat', ...
    'Mihili_CO_FF_2015-06-11.mat', ...
    'Mihili_CO_FF_2015-06-12.mat', ...
    'Mihili_CO_FF_2015-06-15.mat', ...
    'Mihili_CO_FF_2015-06-16.mat', ...
    'Mihili_CO_FF_2015-06-17.mat', ...
    'Mihili_CO_VR_2015-06-23.mat', ...
    'Mihili_CO_VR_2015-06-25.mat', ...
    'Mihili_CO_VR_2015-06-26.mat', ...
    };

files_mrt = { ...
    'MrT_CO_FF_2013-08-19.mat', ...
    'MrT_CO_FF_2013-08-21.mat', ...
    'MrT_CO_FF_2013-08-23.mat', ...
    'MrT_CO_VR_2013-09-03.mat', ...
    'MrT_CO_VR_2013-09-05.mat', ...
    'MrT_CO_VR_2013-09-09.mat', ...
    };



%% get the list of file paths
switch lower(pars.monkey)
    case 'chewie'
                file_list = files_chewie;
    case 'mihili'
                file_list = files_mihili;
    case 'mrt'
        file_list = files_mrt;
    case 'chewie2'
        file_list = files_chewie2;
        pars.monkey  = 'Chewie';
end

clear file_info;
for iFile = 1:length(file_list)
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    
    date = temp{4};
    task = temp{2};
    monkey = temp{1};
    pert = temp{3};
    
    filename = [pars.monkey '_' pars.array '_' task '_' pert '_BL_' datestr(date,'mmddyyyy') '_001.nev'];
    
    %     file_info(iFile).filepath = fullfile(data_root,monkey,'CerebusData',date,filename);
    file_info(iFile).filepath = fullfile(data_root,pars.monkey,filename);
    file_info(iFile).date = date;
    file_info(iFile).task = task;
    file_info(iFile).monkey = monkey;
    file_info(iFile).pert = pert;
end

% sort by date
[~,idx] = sort(cellfun(@(x) datenum(x,'yyyy-mm-dd'),{file_info.date}));
file_info = file_info(idx);


%% pick the sesson comparisons
% get all pairs of sessions
% get all pairs of sessions
sessions                = {file_info.date};
n_sessions              = length(sessions);
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);

diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end


%% process raw data
% Requires access to the raw data
if do_data_processing
    file_data = cell(1,length(file_info));
    for iFile = 1:length(file_info)
        clear data;
        if exist([file_info(iFile).filepath(1:end-4) '.mat'],'file')
            load([file_info(iFile).filepath(1:end-4) '.mat']);
        else
            NEV = openNEV_td(file_info(iFile).filepath,'read','nosave');
        end
        
        % convert NEV to data struct format
        %   cell pars.array where each day is an entry
        %   each day has fields:
        %       .sg (spike guide)
        %       .units(:) (pars.array of units recorded on that day
        %       .units(:).ts (spike times)
        %       .units(:).wf (waveforms)
        %       .units(:).id ([electrode unit] label, duplicate of spike guide
        elecs = unique(NEV.Data.Spikes.Electrode);
        dur = NEV.MetaTags.DataDurationSec;
        
        sg = [];
        c = 0; % counter
        for e = 1:length(elecs)
            if elecs(e) <= 96
                elec_idx = NEV.Data.Spikes.Electrode == elecs(e);
                units = unique(NEV.Data.Spikes.Unit(elec_idx));
                
                units = units(units > 0 & units < 255);
                for u = 1:length(units)
                    
                    unit_idx = NEV.Data.Spikes.Electrode == elecs(e) & NEV.Data.Spikes.Unit == units(u);
                    
                    wf = NEV.Data.Spikes.Waveform(:,unit_idx);
                    ts = NEV.Data.Spikes.TimeStamp(:,unit_idx);
                    id = [elecs(e) units(u)];
                    
                    % calculate some metrics to add in to the data struct, cuz why not
                    mean_wf = mean(wf,1);
                    avg_fr = length(ts)/dur;
                    
                    % package it up
                    c = c+1; % increment counter
                    data.units(c).ts = double(ts)/30000;
                    data.units(c).wf = wf;
                    data.units(c).id = id;
                    data.units(c).mean_wf = mean_wf;
                    data.units(c).avg_fr = avg_fr;
                    sg = [sg; id];
                end
            end
        end
        data.sg = sg;
        data.duration = dur;
        
        file_data{iFile} = data;
    end
    save(fullfile(data_save_dir,[pars.monkey '_' pars.array '_ArrayStabilityData.mat']),'file_data','file_info','-v7.3');
end

%% do statistical test
% TAKES FOREVER TO RUN!!!!
if do_tracking
    load(fullfile(data_save_dir,[pars.monkey '_' pars.array '_ArrayStabilityData.mat']));
    
    [COMPS, ts_ISI, D_wave, lda_proj] = KS_p({'isi','wf'},file_data,0.95);
    
    % save results
    save(fullfile(out_dir,[pars.monkey '_' pars.array '_ArrayStabilityResults.mat']),'COMPS','file_info');
    curr_files = file_info;
else
    curr_files = file_info;
    load(fullfile(out_dir,[pars.monkey '_' pars.array '_ArrayStabilityResults.mat']),'COMPS','file_info');
    
end


%% find the percent of stable cells for all comparisons

perc_stable = zeros(n_comb_sessions,1);
for c = 1:n_comb_sessions
    % use file_info to match the dates
    old_dates = {file_info.date};
    session1 = find(strcmpi(old_dates,curr_files(comb_sessions(c,1)).date));
    session2 = find(strcmpi(old_dates,curr_files(comb_sessions(c,2)).date));;
    
    % how many neurons from day 'i' ("first day") are matched on Day j
    % (percent is found as divided by number of neurons on day i)
    perc_stable(c) = 100*sum(COMPS{session2}.chan(:,session1) ~= 0)/size(COMPS{session1}.chan,1);
end


%% package up the results
if save_results
    results = struct( ...
        'perc_stable',perc_stable, ...
        'file_info',file_info, ...
        'diff_days',diff_days, ...
        'comb_sessions',comb_sessions, ...
        'pars',pars);
    
    fn = fullfile(out_dir,['SortedNeuronStability_' pars.monkey '_' pars.array '.mat']);
    save(fn,'results');
end



%% plot the stuff
figure; hold all;
plot(diff_days,perc_stable,'ko','LineWidth',2);

% fit a line
[b,~,~,~,s] = regress(perc_stable,[ones(size(diff_days))' diff_days']);
r = s(1);
p = s(3);

V = axis;
plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(diff_days)+1],'YLim',[0 100]);
xlabel('Days between sessions');
ylabel('% of Day 1 neurons that match');

title([pars.monkey ' - ' pars.array]);

if save_figs
    fn = fullfile(out_dir,['SortedNeuronStability_' pars.monkey '_' pars.array]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end

