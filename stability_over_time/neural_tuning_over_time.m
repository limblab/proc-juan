clear;
clc;
close all;

data_dir = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/';
fig_dir = '/Users/mattperich/Dropbox/Research/Papers/2018 - Stability latent activity/RawFigs/';
save_figs = true;

monkey = 'Han';
array = 'S1';
which_tuning_params = [1,2,3]; %1 = mean firing rate, 2 = modulation depth, 3 = preferred direction

min_r = 0.5;

% simple tuning script
pars.bad_neuron_params.min_fr           = 0.1;
% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;

switch lower(array)
    case 'm1'
        pars.tuning_window = {'idx_go_cue',{'idx_go_cue',40}};
    case 'pmd'
        pars.tuning_window = {{'idx_go_cue',-40},'idx_go_cue'};
    case 's1'
        pars.tuning_window = {'idx_goCueTime',{'idx_goCueTime',40}};
end

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

files_han = { ...
    'Han_CO_BB_2017-10-24.mat', ...
    'Han_CO_BB_2017-10-30.mat', ...
    'Han_CO_BB_2017-10-31.mat', ...
    'Han_CO_BB_2017-11-03.mat', ...
    'Han_CO_BB_2017-11-16.mat', ...
    'Han_CO_BB_2017-11-20.mat', ...
    'Han_CO_BB_2017-11-21.mat', ...
    'Han_CO_BB_2017-11-22.mat', ...
    'Han_CO_BB_2017-11-27.mat', ...
    'Han_CO_BB_2017-11-28.mat', ...
    'Han_CO_BB_2017-11-29.mat', ...
    'Han_CO_BB_2017-12-01.mat', ...
    'Han_CO_BB_2017-12-04.mat', ...
    'Han_CO_BB_2017-12-07.mat', ...
    };


%% get the list of file paths
switch lower(monkey)
    case 'chewie'
        file_list = files_chewie;
    case 'mihili'
        file_list = files_mihili;
    case 'han'
        if ~strcmpi(array,'s1')
            error('Han only has an S1 array.');
        end
        file_list = files_han;
end

clear file_info;
for iFile = 1:length(file_list)
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    
    date = temp{4};
    task = temp{2};
    monkey = temp{1};
    pert = temp{3};
    
    %     file_info(iFile).filepath = fullfile(data_root,monkey,'CerebusData',date,filename);
    file_info(iFile).date = date;
    file_info(iFile).task = task;
    file_info(iFile).monkey = monkey;
    file_info(iFile).pert = pert;
end

% sort by date
[~,idx] = sort(cellfun(@(x) datenum(x,'yyyy-mm-dd'),{file_info.date}));
file_info = file_info(idx);

if strcmpi(monkey,'han')
    load(fullfile(data_dir,'han_tds.mat'));
    
    load_td = [];
    dates  = unique({trial_data.date});
    for iDate = 1:length(dates)
        [~,td] = getTDidx(trial_data,'result','R','date',dates{iDate});
        td = td(~isnan([td.idx_goCueTime]));
        td = stripSpikeSorting(td);
        td = removeBadNeurons(td,pars.bad_neuron_params);
        td = trimTD(td,pars.tuning_window{1},pars.tuning_window{2});
        
        load_td = [load_td, td];
    end
    trial_data = load_td; clear load_td td;
    
else
    filepaths = cell(size(file_list));
    for iFile = 1:length(file_list)
        filepaths{iFile} = fullfile(data_dir,file_list{iFile});
    end
    load_fn_call = { ...
        {@getTDidx,'result','R'}, ...
        @stripSpikeSorting, ...
        {@removeBadNeurons,pars.bad_neuron_params}, ...
        {@trimTD,pars.tuning_window{1},pars.tuning_window{2}}, ...
        };
    trial_data = loadTDfiles(filepaths, load_fn_call{:});
end

%% do the tuning
file_sgs = cell(length(file_list),1);
file_tuning = cell(length(file_list),3);
for iFile = 1:length(file_list)
    % get the data for this file
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    the_date = temp{4};
    
    [~,td] = getTDidx(trial_data,'date',datestr(the_date,'mm-dd-yyyy'));
    
    % get average firing rate in window
    fr = cell2mat(cellfun(@(x) sum(x,1)/(0.01*size(x,1)),{td.([array '_spikes'])},'uni',0)');
    % get target directions
    dir = [td.target_direction];
    % fit tuning curves
    [tc,cb,r] = regressTuningCurves(fr,dir',{'bootstrap',100,0.95});
    % store it for later
    file_tuning{iFile,1} = tc;
    file_tuning{iFile,2} = cb;
    file_tuning{iFile,3} = r;
    file_sgs{iFile} = td(1).([array '_unit_guide']);
    
    
end


%% plot the stuff

for tp = 1:length(which_tuning_params)
    
    switch which_tuning_params(tp)
        case 1
            y_label = '|Change in Mean Firing| (Hz)';
            file_addendum = 'MeanFiring';
        case 2
            y_label = '|Change in Modulation Depth| (Hz)';
            file_addendum = 'ModDepth';
        case 3
            y_label = '|Change in Preferred Direction| (Deg)';
            file_addendum = 'PD';
    end
    
    
    figure;
    hold all;
    
    [time_diff, tune_diff] = deal([]);
    for i = 1:size(file_tuning,1)-1
        for j = i+1:size(file_tuning,1)
            % make sure you study the same elctrodes on each day
            [elecs,idx1,idx2] = intersect(file_sgs{i}(:,1),file_sgs{j}(:,1));
            
            % make sure the electrodes are tuned on each day
            r1 = file_tuning{i,3}(idx1);
            r2 = file_tuning{j,3}(idx2);
            
            idx1 = idx1( r1 > min_r & r2 > min_r );
            idx2 = idx2( r1 > min_r & r2 > min_r );
            
            if which_tuning_params(tp) == 3
                d = angleDiff(file_tuning{i,1}(idx1,3),file_tuning{j,1}(idx2,3),true,false);
                m = 180/pi*circular_mean(d);
                s = 180/pi*circular_std(d)/sqrt(length(d));
            else
                d = abs(file_tuning{j,1}(idx2,which_tuning_params(tp)) - file_tuning{i,1}(idx1,which_tuning_params(tp)));
                m = mean(d);
                s = std(d)/sqrt(length(d));
            end
            x = datenum(file_info(j).date,'yyyy-mm-dd') - datenum(file_info(i).date,'yyyy-mm-dd');
            
            plot(x,m,'ko','LineWidth',2);
            plot([x, x],[m-s,m+s],'k-','LineWidth',2)
            
            time_diff = [time_diff, x];
            tune_diff = [tune_diff, m];
        end
    end
    
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(time_diff)+1]);
    
    % fit a line
    [b,~,~,~,s] = regress(tune_diff',[ones(size(time_diff))' time_diff']);
    r = s(1);
    p = s(3);
    
    V = axis;
    plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
    text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);
    
    
    
    title([monkey ' - ' array]);
    xlabel('Days between sessions')
    ylabel(y_label);
    
    
    if save_figs
        fn = fullfile(fig_dir,['NeuralTuningStability_' file_addendum '_' monkey '_' array]);
        savefig(gcf,[fn '.fig']);
        saveas(gcf,fn,'png');
        saveas(gcf,fn,'pdf');
    end
    
end
