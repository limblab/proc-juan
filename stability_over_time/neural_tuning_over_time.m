clear;
clc;
close all;

data_dir = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/';
out_dir = '/Users/mattperich/Dropbox/Research/Papers/Juan and Matt - Stability latent activity/Results/Tuning stability/';
save_figs = true;
save_results = true;
do_the_boxplots = false;

pars.monkey = 'Mihili';
pars.array = 'PMd';

% will plot all of these
which_tuning_params = [1,2,3]; %1 = mean firing rate, 2 = modulation depth, 3 = preferred direction
pars.tuning_param_labels = {'Mean Firing Rate','Modulation Depth','Preferred Direction'};

% use this to throw away some super-badly-tuned cells just to denoise a bits
%   note: currently, the lower 95% C.I. of R2 must be above this value
pars.min_r2 = 0;

% I just use this match Juan's tuning window calculation... in reality I
% don't bother ACTUALLY downsampling
pars.n_bins_downs = 3;

% simple tuning script
pars.bad_neuron_params.min_fr           = 0.1;
% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;

% I did 1000 for the saved figures, but 100 is good enough in general
pars.tuning_conf_params = {'bootstrap',100,0.95};

switch lower(pars.monkey)
    case {'han','chips'}
        pars.idx_start          = {'idx_goCueTime',round(-2*5/3 * pars.n_bins_downs)};
        pars.idx_end            = {'idx_goCueTime',round(9*5/3 * pars.n_bins_downs)};
    case {'chewie','mihili','jaco','chewie2','mrt'}
        switch lower(pars.array)
            case 'm1'
                pars.idx_start  = {'idx_movement_on',-2 * pars.n_bins_downs};
                pars.idx_end    = {'idx_movement_on',13 * pars.n_bins_downs};
            case 'pmd'
                pars.idx_start  = {'idx_go_cue',-13 * pars.n_bins_downs};
                pars.idx_end    = {'idx_go_cue',2 * pars.n_bins_downs};
        end
end

monkey_lists;


%% get the list of file paths
switch lower(pars.monkey)
    case 'chewie'
        file_list = files_chewie;
    case 'chewie2'
        file_list = files_chewie2;
    case 'mihili'
        file_list = files_mihili;
    case 'mrt'
        file_list = files_mrt;
    case 'jaco'
        file_list = files_jaco;
    case 'chips'
        if ~strcmpi(pars.array,'s1')
            error('Han only has an S1 array.');
        end
        file_list = files_chips;
    case 'han'
        if ~strcmpi(pars.array,'s1')
            error('Han only has an S1 array.');
        end
        file_list = files_han;
end

clear file_info;
for iFile = 1:length(file_list)
    [~,temp,~] = fileparts(file_list{iFile});
    temp = strsplit(temp,'_');
    
    % there are all kinds of different date formats and missing information
    % and shit, so I needed to add some special cases
    if strcmpi(pars.monkey,'chips')
        date = datestr(datenum(temp{2},'yyyymmdd'),'mm-dd-yyyy');
        task = 'CO';
        monkey = 'Chips';
        pert = '';
    else
        date = datestr(datenum(temp{4},'yyyy-mm-dd'),'mm-dd-yyyy');
        task = temp{2};
        monkey = temp{1};
        pert = temp{3};
    end
    
    %     file_info(iFile).filepath = fullfile(data_dir,monkey,'CerebusData',date,filename);
    file_info(iFile).date = date;
    file_info(iFile).task = task;
    file_info(iFile).monkey = monkey;
    file_info(iFile).pert = pert;
end

% sort by date
[~,idx] = sort(cellfun(@(x) datenum(x,'mm-dd-yyyy'),{file_info.date}));
file_info = file_info(idx);

%% LOAD THE DATA
if strcmpi(pars.monkey,'han')
    load(fullfile(data_dir,'han_tds.mat'));
    
    load_td = [];
    dates  = unique({trial_data.date});
    for iFile = 1:length(dates)
        [~,td] = getTDidx(trial_data,'result','R','date',dates{iFile});
        td = td(~isnan([td.idx_goCueTime]));
        td = stripSpikeSorting(td);
        td = removeBadTrials(td);
        td = removeBadNeurons(td,pars.bad_neuron_params);
        td = trimTD(td,pars.idx_start,pars.idx_end);
        
        load_td = [load_td, td];
    end
    trial_data = load_td; clear load_td td;
    
elseif strcmpi(pars.monkey,'chips')
    
    load_td = [];
    for iFile = 1:length(file_list)
        load(fullfile(data_dir,file_list{iFile}));
        [~,td] = getTDidx(trial_data,'result','R');
        td = td(~isnan([td.idx_goCueTime]));
        td = stripSpikeSorting(td);
        td = removeBadTrials(td);
        td = removeBadNeurons(td,pars.bad_neuron_params);
        td = trimTD(td,pars.idx_start,pars.idx_end);
        
        load_td = [load_td, td];
    end
    trial_data = load_td; clear load_td td;
    
else % chewie and mihili both fall in this one
    filepaths = cell(size(file_list));
    for iFile = 1:length(file_list)
        filepaths{iFile} = fullfile(data_dir,file_list{iFile});
    end
    load_fn_call = { ...
        {@getTDidx,'result','R','epoch','BL'}, ...
        @stripSpikeSorting, ...
        @removeBadTrials, ...
        {@removeBadNeurons,pars.bad_neuron_params}, ...
        {@trimTD,pars.idx_start,pars.idx_end}, ...
        };
    trial_data = loadTDfiles(filepaths, load_fn_call{:});
end


%% pick the sesson comparisons
% get all pairs of sessions
% get all pairs of sessions
sessions                = unique({trial_data.date});
n_sessions              = length(sessions);
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);

diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end


%% do the tuning
file_sgs = cell(length(file_list),1);
file_tuning = cell(length(file_list),3);
for iFile = 1:length(file_list)
    % get the data for this file
    the_date = file_info(iFile).date;
    
    [~,td] = getTDidx(trial_data,'date',the_date);
    
    % get average firing rate in window
    fr = cell2mat(cellfun(@(x) sum(x,1)/(0.01*size(x,1)),{td.([pars.array '_spikes'])},'uni',0)');
    % get target directions
    dir = [td.target_direction];
    % fit tuning curves
    [tc,cb,r] = regressTuningCurves(fr,dir',pars.tuning_conf_params);
    
    % get confidence bounds for bootstrapped r2 vals
    r = prctile(r,[2.5 97.5],2);
    
    % store it for later
    file_tuning{iFile,1} = tc;
    file_tuning{iFile,2} = cb;
    file_tuning{iFile,3} = r;
    file_sgs{iFile} = td(1).([pars.array '_unit_guide']);
    
    
end


%% do the comparison, save the results
[tuneDiff_m, tuneDiff_sd] = deal(zeros(n_comb_sessions,length(which_tuning_params)));
tuneDiff_all = cell(n_comb_sessions,length(which_tuning_params));
for tp = 1:length(which_tuning_params)
    
    for c = 1:n_comb_sessions
        % make sure you study the same elctrodes on each day
        [elecs,idx1,idx2] = intersect(file_sgs{comb_sessions(c,1)}(:,1),file_sgs{comb_sessions(c,2)}(:,1));
        
        % make sure the electrodes are tuned on each day
        r1 = file_tuning{comb_sessions(c,1),3}(idx1,1);
        r2 = file_tuning{comb_sessions(c,2),3}(idx2,1);
        
        idx1 = idx1( r1 > pars.min_r2 & r2 > pars.min_r2 );
        idx2 = idx2( r1 > pars.min_r2 & r2 > pars.min_r2 );
        
        if which_tuning_params(tp) == 3 % PD, do some angular stuff
            d = angleDiff(file_tuning{comb_sessions(c,1),1}(idx1,3),file_tuning{comb_sessions(c,2),1}(idx2,3),true,false);
            m = 180/pi*circular_mean(d);
            s = 180/pi*circular_std(d)/sqrt(length(d));
            d = 180/pi*d; % for tuneDiff_all, convert it to degrees
        else % mean firing and modulation depth
            d = abs(file_tuning{comb_sessions(c,2),1}(idx2,which_tuning_params(tp)) - file_tuning{comb_sessions(c,1),1}(idx1,which_tuning_params(tp)));
            m = mean(d);
            s = std(d)/sqrt(length(d));
        end
        
        tuneDiff_m(c,tp) = m;
        tuneDiff_sd(c,tp) = s;
        tuneDiff_all{c,tp} = d;
    end
    
    
end

%% package up dem results
if save_results
    results = struct( ...
        'tuneDiff_m',tuneDiff_m, ...
        'tuneDiff_sd',tuneDiff_sd, ...
        'tuneDiff_all',{tuneDiff_all}, ...
        'file_info',file_info, ...
        'file_tuning',{file_tuning}, ...
        'file_sgs',{file_sgs}, ...
        'diff_days',diff_days, ...
        'comb_sessions',comb_sessions, ...
        'pars',pars);
    
    fn = fullfile(out_dir,['NeuralTuningStability_' pars.monkey '_' pars.array '.mat']);
    save(fn,'results');
end



%% now plot stuff

figure('Position',[100 100 800 400]);
for tp = 1:length(which_tuning_params)
    
    subplot(1,length(which_tuning_params),tp);
    hold all;
    
    [d,g] = deal([]);
    for c = 1:n_comb_sessions
        m = tuneDiff_m(c,tp);
        s = tuneDiff_sd(c,tp);
        
        x = datenum(file_info(comb_sessions(c,2)).date,'mm-dd-yyyy') - datenum(file_info(comb_sessions(c,1)).date,'mm-dd-yyyy');
        
        plot(x,m,'ko','LineWidth',2);
        plot([x, x],[m-s,m+s],'k-','LineWidth',2)
        
    end
    
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(diff_days)+1]);
    
    % fit a line
    [b,~,~,~,s] = regress(tuneDiff_m(:,tp),[ones(size(diff_days))' diff_days']);
    r = s(1);
    p = s(3);
    
    V = axis;
    plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
    text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);
    
    title([pars.monkey ' - ' pars.array]);
    xlabel('Days between sessions')
    ylabel(['|Change in ' pars.tuning_param_labels{tp} '|']);
    
end

if save_figs
    fn = fullfile(out_dir,['NeuralTuningStability_' pars.monkey '_' pars.array]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end

%% this version does the boxplot
if do_the_boxplots
    figure('Position',[100 100 1200 400]);
    for tp = 1:length(which_tuning_params)
        
            subplot(1,length(which_tuning_params),tp);
        hold all;
        
        [d,g] = deal([]);
        for c = 1:n_comb_sessions
            m = tuneDiff_m(c,tp);
            s = tuneDiff_sd(c,tp);
            
            x = datenum(file_info(comb_sessions(c,2)).date,'mm-dd-yyyy') - datenum(file_info(comb_sessions(c,1)).date,'mm-dd-yyyy');
            
            d = [d; tuneDiff_all{c,tp}];
            g = [g; x*ones(size(tuneDiff_all{c,tp}))];
            
            %         plot(x,m,'ko','LineWidth',2);
            %         plot([x, x],[m-s,m+s],'k-','LineWidth',2)
            
        end
        boxplot(d,g);

        set(gca,'Box','off','TickDir','out','FontSize',14);
        
        % fit a line
        [b,~,~,~,s] = regress(tuneDiff_m(:,tp),[ones(size(diff_days))' diff_days']);
        r = s(1);
        p = s(3);
        
        V = axis;
        plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);
        text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);
        
        title([pars.monkey ' - ' pars.array]);
        xlabel('Days between sessions')
        ylabel(['|Change in ' pars.tuning_param_labels{tp} '|']);
        
        if save_figs
            fn = fullfile(out_dir,['NeuralTuningStability_' pars.monkey '_' pars.array '_Boxplot']);
            savefig(gcf,[fn '.fig']);
            saveas(gcf,fn,'png');
            saveas(gcf,fn,'pdf');
        end
    end
end