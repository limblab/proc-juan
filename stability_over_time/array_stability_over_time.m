clear;
clc;
close all;

data_root = '/Volumes/MattData/';
data_save_dir = '/Users/mattperich/Dropbox/Research/Data/StabilityData';
out_dir = '/Users/mattperich/Dropbox/Research/Papers/Juan and Matt - Stability latent activity/Results/Neuron stability/';

% data_root = 'C:\Users\Matt\Desktop\LimblabData\';
% data_save_dir = 'C:\Users\Matt\Desktop\LimblabData\';
% out_dir = 'C:\Users\Matt\Desktop\LimblabData\';

do_data_processing = false; % can only do this if you have raw data
do_tracking        = false; % can only do this if you have completed data processing
save_results       = true; % save a .mat file with the results
save_figs          = true;

criteria = {'isi','wave'};
p_val = 0.1;
do_norm = true;


pars.monkey = 'Chewie2';
pars.array = 'M1';
pars.spiking_inputs{1} = [pars.array '_spikes'];

monkey_file_lists;


%% get the list of file paths
switch lower(pars.monkey)
    case 'chewie'
        file_list = files_chewie;
        use_monkey_name = 'Chewie';
    case 'mihili'
        file_list = files_mihili;
        use_monkey_name = 'Mihili';
    case 'mrt'
        file_list = files_mrt;
        use_monkey_name = 'MrT';
    case 'chewie2'
        file_list = files_chewie2;
        use_monkey_name = 'Chewie';
end

clear file_info;
for iFile = 1:length(file_list)
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    
    date = temp{4};
    task = temp{2};
    monkey = temp{1};
    pert = temp{3};
    
    filename = [use_monkey_name '_' pars.array '_' task '_' pert '_BL_' datestr(date,'mmddyyyy') '_001.nev'];
    
    %     file_info(iFile).filepath = fullfile(data_root,monkey,'CerebusData',date,filename);
    file_info(iFile).filepath = fullfile(data_root,use_monkey_name,filename);
    file_info(iFile).date = date;
    file_info(iFile).task = task;
    file_info(iFile).monkey = monkey;
    file_info(iFile).pert = pert;
    
end

% sort by date
[~,idx] = sort(cellfun(@(x) datenum(x,'yyyy-mm-dd'),{file_info.date}));
file_info = file_info(idx);

% get time between days
date_diff = cellfun(@datenum,{file_info.date});
date_diff = date_diff - min(date_diff);


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
    load(fullfile(data_save_dir,'StabilityData',[pars.monkey '_' pars.array '_ArrayStabilityData.mat']));
    
    [COMPS, ts_ISI, D_wave, lda_proj] = KS_p({'isi','wf'},file_data,(1-p_val));
    
    % save results
    save(fullfile(out_dir,[pars.monkey '_' pars.array '_ArrayStabilityResults.mat']),'COMPS','file_info');
    curr_files = file_info;
else
    curr_files = file_info;
    load(fullfile(out_dir,[pars.monkey '_' pars.array '_ArrayStabilityResults.mat']),'COMPS','file_info');
    
end


%% calculate the matches

% compare units
num_days =  length(COMPS);

all_match = cell(1,num_days);
for i = 1:num_days % Loop through days
    num_cells = size(COMPS{i}.chan,1);
    
    match_days = zeros(num_cells,num_days);
    for j = 1:num_cells % find ID of unit in day i
        for k = find(1:num_days ~= i) % Look at other days
            % Find units on the same channel
            all_p_isi = COMPS{i}.p_isi{j,k};
            all_p_wave = COMPS{i}.p_wave{j,k};
            
            match_temp = zeros(1,size(all_p_isi,1));
            for l = 1:size(all_p_isi,1) % For all units on the same electrode
                p_isi = all_p_isi(l,2);
                p_wave = all_p_wave(l,2);
                
                temp = 1;
                if ismember('isi',criteria)
                    temp = temp*p_isi;
                end
                if ismember('wave',criteria)
                    temp = temp*p_wave;
                end
                
                if temp < (p_val)^length(criteria)
                    match_temp(l) = 1;
                end
                
            end
            
            if sum(match_temp) > 1
                %                 error('Multiple matches OH SHIT');
            end
            match_days(j,k) = sum(match_temp);
        end
    end
    
    all_match{i} = match_days;
end


% find the percent of stable cells for all comparisons

[perc_stable,num_cells] = deal(zeros(n_comb_sessions,1));
for c = 1:n_comb_sessions
    % use file_info to match the dates
    old_dates = {file_info.date};
    session1 = find(strcmpi(old_dates,curr_files(comb_sessions(c,1)).date));
    session2 = find(strcmpi(old_dates,curr_files(comb_sessions(c,2)).date));

%     if datenum(curr_files(comb_sessions(c,1)).date,'yyyy-mm-dd') == datenum('20161005','yyyymmdd') && ...
%             datenum(curr_files(comb_sessions(c,2)).date,'yyyy-mm-dd') == datenum('20161021','yyyymmdd')
%         keepit = c; 
%     end
    
    num_cells(c) = min( [size(COMPS{session1}.chan,1), size(COMPS{session2}.chan,1) ] );
    
    % how many neurons from day 'i' ("first day") are matched on Day j
    % (percent is found as divided by number of neurons on day i)
    if do_norm
        perc_stable(c) = 100*sum(all_match{session2}(:,session1) ~= 0)/num_cells(c);
    else
        perc_stable(c) = sum(all_match{session2}(:,session1) ~= 0);
    end
end


%% package up the results
if save_results
    results = struct( ...
        'perc_stable',perc_stable, ...
        'file_info',file_info, ...
        'curr_files',curr_files, ...
        'diff_days',diff_days, ...
        'comb_sessions',comb_sessions, ...
        'COMPS',COMPS, ...
        'pars',pars);
    
    fn = fullfile(out_dir,['SortedNeuronStability_' pars.monkey '_' pars.array '.mat']);
    save(fn,'results');
end



%% plot the stuff
figure; hold all;
plot(diff_days,perc_stable,'.','markersize',32,'color',[0.3 0.3 0.3]);

% fit a line
[b,~,~,~,s] = regress(perc_stable,[ones(size(diff_days))' diff_days']);
r = s(1);
p = s(3);

V = axis;
plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(diff_days)+1],'YLim',[0 100]);
xlabel('Days between sessions');
if do_norm
    ylabel('% of neurons that match');
else
    ylabel('# of neurons that match');
end

title([pars.monkey ' - ' pars.array]);

if save_figs
    if do_norm
        fn = fullfile(out_dir,[pars.monkey '_' pars.array '_SortedNeuronStability_Percent' ]);
    else
        fn = fullfile(out_dir,[pars.monkey '_' pars.array '_SortedNeuronStability' ]);
    end
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end


%%
figure; hold all;
old_dates = {file_info.date};
for i = 1:length(curr_files)
    session1 = find(strcmpi(old_dates,curr_files(i).date));
    plot(date_diff(i),size(COMPS{session1}.chan,1),'.','markersize',32,'color',[0.3 0.3 0.3]);
end
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Days since first recording');
ylabel('Number of sorted neurons');

title([pars.monkey ' - ' pars.array]);

if save_figs
    fn = fullfile(out_dir,[pars.monkey '_' pars.array '_SortedNeuronCounts' ]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end


