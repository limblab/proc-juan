dclear;
clc;
close all;

redo_data = false;
redo_results = false;

data_root = '/Volumes/MattData/';
save_dir = '/Users/mattperich/Dropbox/Research/Data/StabilityData';

monkey = 'Chewie';
array = 'M1';

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
    'Chewie_CO_FF_2016-10-13.mat' ...
    };

files_mihili = { ...
    'Mihili_CO_FF_2014-02-03.mat', ...
    'Mihili_CO_FF_2014-02-17.mat', ...
    'Mihili_CO_FF_2014-02-18.mat', ...
    'Mihili_CO_VR_2014-03-03.mat', ...
    'Mihili_CO_VR_2014-03-04.mat', ...
    'Mihili_CO_VR_2014-03-06.mat', ...
    'Mihili_CO_FF_2014-03-07.mat' ...
    };

% took these parameters from Juan
bad_neuron_params.min_fr           = 0.1;
bad_neuron_params.shunt_check_yn   = false;


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
    
    file_info(iFile).filepath = fullfile(data_root,monkey,'CerebusData',date,filename);
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
        NEV = openNEV_td(file_info(iFile).filepath,'read','nosave');
        
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
            elec_idx = NEV.Data.Spikes.Electrode == elecs(e);
            units = unique(NEV.Data.Spikes.Unit(elec_idx));
            
            units = units(units > 0 & units < 255);
            for u = 1:length(units)
                c = c+1; % increment counter
                
                unit_idx = NEV.Data.Spikes.Electrode == elecs(e) & NEV.Data.Spikes.Unit == units(u);
                
                wf = NEV.Data.Spikes.Waveform(:,unit_idx);
                ts = NEV.Data.Spikes.TimeStamp(:,unit_idx);
                id = [elecs(e) units(u)];
                
                % calculate some metrics to add in to the data struct, cuz why not
                mean_wf = mean(wf,1);
                avg_fr = length(ts)/dur;
                
                % package it up
                data.units(c).ts = double(ts)/30000;
                data.units(c).wf = wf;
                data.units(c).id = id;
                data.units(c).mean_wf = mean_wf;
                data.units(c).avg_fr = avg_fr;
                sg = [sg; id];
            end
        end
        data.sg = sg;
        data.duration = dur;
        
        file_data{iFile} = data;
    end
    save(fullfile(save_dir,[monkey '_' array '_ArrayStabilityData.mat']),'file_data','file_info');
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

temp = zeros(1,length(COMPS));
for i = 1:length(COMPS)
    temp(i) = sum(COMPS{i}.chan(:,1) ~= 0)/size(COMPS{1}.chan,1);
    the_dates(i) = datenum(file_info(i).date,'yyyy-mm-dd');
end

figure;
plot(the_dates - the_dates(1),100*temp,'LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0 100]);
xlabel('Days since first session');
ylabel('% of Day 1 neurons that match');

title([monkey ' - ' array]);
%%
% count = 0;
% for iDay = 1:size(goodDates,1)
%
%     comp = COMPS{iDay}.chan;
%     for unit = 1:size(comp,1)
%         % make sure this unit hasn't been used in the past
%         %   loop along all previous days and make sure it hasn't matched
%         newCheck = zeros(1,iDay-1);
%         for j = 1:iDay-1
%             % get matching info and index of unit ID
%             comp2 = tracking.(useArray){j}.chan;
%             idx = find(comp2(:,j)==comp(unit,iDay));
%
%             % if the cell has already matched with this or a later one
%             if ~isempty(idx)
%                 newCheck(j) = any(comp2(idx,iDay:end) > 0);
%             end
%         end
%
%         % if no previous day matched with this cell
%         if ~any(newCheck)
%             count = count+1;
%
%             e = floor(comp(unit,iDay));
%             u = int32(10*rem(comp(unit,iDay),e));
%
%             ind = sg(:,1)==e & sg(:,2) == u;
%
%             for j = iDay:size(comp,2)
%                 if comp(unit,j)
%
%                     e = floor(comp(unit,j));
%                     u = int32(10*rem(comp(unit,j),e));
%
%                     idx = sg(:,1)==e & sg(:,2)==u;
%                 end
%             end
%         end
%
%     end
% end
%
