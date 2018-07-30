clear;
clc;
close all;



% data_root = '/Volumes/MattData/';
data_root = 'C:\Users\Matt\Desktop\LimblabData';
% save_dir = '/Users/mattperich/Dropbox/Research/Data/StabilityData';
save_dir = 'C:\Users\Matt\Desktop\LimblabData\StabilityData';%'C:\Users\Matt\Dropbox\Research\Data\StabilityData';
fig_dir = '/Users/mattperich/Dropbox/Research/Papers/2018 - Stability latent activity/Figs/raw/';

redo_data = false;
redo_results = true;
save_figs = false;

monkey = 'Chewie';
array = 'PMd';


files_chewie = { ...
    'Chewie_CO_VR_2016-09-12.mat', ...
    'Chewie_CO_VR_2016-09-14.mat', ...
    'Chewie_CO_FF_2016-09-15.mat', ...
    'Chewie_CO_FF_2016-09-19.mat', ...
    'Chewie_CO_FF_2016-09-21.mat', ...
    'Chewie_CO_FF_2016-10-05.mat', ...
    'Chewie_CO_VR_2016-10-06.mat', ...
    'Chewie_CO_FF_2016-10-07.mat', ...
    'Chewie_CO_FF_2016-10-11.mat', ...
    'Chewie_CO_FF_2016-10-13.mat' ...% % add these below here to TD maybe
    'Chewie_CO_VR_2016-09-09.mat', ...
    'Chewie_CO_FF_2016-09-23.mat', ...
    'Chewie_CO_VR_2016-09-29.mat', ...
    'Chewie_CO_CS_2016-10-14.mat', ... % this one was sorted without the Ricardo split
    'Chewie_CO_CS_2016-10-21.mat', ...
    };

files_mihili = { ...
    'Mihili_CO_FF_2014-02-03.mat', ...
    'Mihili_CO_FF_2014-02-17.mat', ...
    'Mihili_CO_FF_2014-02-18.mat', ...
    'Mihili_CO_VR_2014-03-03.mat', ...
    'Mihili_CO_VR_2014-03-04.mat', ...
    'Mihili_CO_VR_2014-03-06.mat', ...
    'Mihili_CO_FF_2014-03-07.mat', ...
    'Mihili_CO_FF_2015-06-10.mat', ...
    'Mihili_CO_FF_2015-06-11.mat', ...
    'Mihili_CO_FF_2015-06-15.mat', ...
    'Mihili_CO_FF_2015-06-16.mat', ...
    'Mihili_CO_FF_2015-06-17.mat', ...
    'Mihili_CO_VR_2015-06-23.mat', ...
    'Mihili_CO_VR_2015-06-25.mat', ...
    'Mihili_CO_VR_2015-06-26.mat', ...
    };


%% get the list of file paths
switch lower(monkey)
    case 'chewie'
        file_list = files_chewie;
    case 'mihili'
        file_list = files_mihili;
end

clear file_info;
for iFile = 1:length(file_list)
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    
    date = temp{4};
    task = temp{2};
    monkey = temp{1};
    pert = temp{3};
    
    filename = [monkey '_' array '_' task '_' pert '_BL_' datestr(date,'mmddyyyy') '_001.nev'];
    
    %     file_info(iFile).filepath = fullfile(data_root,monkey,'CerebusData',date,filename);
    file_info(iFile).filepath = fullfile(data_root,filename);
    file_info(iFile).date = date;
    file_info(iFile).task = task;
    file_info(iFile).monkey = monkey;
    file_info(iFile).pert = pert;
end


%% process raw data
if ~exist(fullfile(save_dir,[monkey '_' array '_ArrayStabilityData.mat']),'file') || redo_data
    file_data = cell(1,length(file_info));
    for iFile = 1:length(file_info)
        clear data;
        if exist([file_info(iFile).filepath(1:end-4) '.mat'],'file')
            load([file_info(iFile).filepath(1:end-4) '.mat']);
        else
            NEV = openNEV_td(file_info(iFile).filepath,'read','nosave');
        end
        
        % convert NEV to data struct format
        %   cell array where each day is an entry
        %   each day has fields:
        %       .sg (spike guide)
        %       .units(:) (array of units recorded on that day
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
    save(fullfile(save_dir,[monkey '_' array '_ArrayStabilityData.mat']),'file_data','file_info','-v7.3');
end

%% do statistical test
if ~exist(fullfile(save_dir,[monkey '_' array '_ArrayStabilityResults.mat']),'file') || redo_results
    load(fullfile(save_dir,[monkey '_' array '_ArrayStabilityData.mat']));
    
    [COMPS, ts_ISI, D_wave, lda_proj] = KS_p({'isi','wf'},file_data,0.95);
    
    % save results
    save(fullfile(save_dir,[monkey '_' array '_ArrayStabilityResults.mat']),'COMPS','ts_ISI','D_wave','lda_proj');
else
    load(fullfile(save_dir,[monkey '_' array '_ArrayStabilityResults.mat']),'COMPS','ts_ISI','D_wave','lda_proj');
    
end


%%

perc_stable = [];
time_diff = [];
for i = 1:length(COMPS)-1
    for j = i+1:length(COMPS)
        % how many neurons from day 'i' ("first day") are matched on Day j
        % (percent is found as divided by number of neurons on day i)
        perc_stable = [perc_stable, 100*sum(COMPS{j}.chan(:,i) ~= 0)/size(COMPS{i}.chan,1)];
        time_diff = [time_diff, datenum(file_info(j).date,'yyyy-mm-dd') - datenum(file_info(i).date,'yyyy-mm-dd')];
    end
end

figure; hold all;
plot(time_diff,perc_stable,'ko','LineWidth',2);

% fit a line
[b,~,~,~,s] = regress(perc_stable',[ones(size(time_diff))' time_diff']);
r = s(1);
p = s(3);

V = axis;
plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(time_diff)+1],'YLim',[0 100]);
xlabel('Days between sessions');
ylabel('% of Day 1 neurons that match');

title([monkey ' - ' array]);

if save_figs
    fn = fullfile(fig_dir,['SortedNeuronStability' '_' monkey '_' array]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end

