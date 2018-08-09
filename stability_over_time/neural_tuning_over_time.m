clear;
clc;
close all;

[save_dir, data_dir] = get_computer_paths();

save_figs = true;
save_results = true;
do_the_boxplots = false;

pars.monkey = 'Chewie';
pars.array = 'M1';
pars.spiking_inputs{1} = [pars.array '_spikes'];

% will plot all of these
which_tuning_params = [1,2,3]; %1 = mean firing rate, 2 = modulation depth, 3 = preferred direction
pars.tuning_param_labels = {'Mean Firing Rate','Modulation Depth','Preferred Direction'};

% use this to throw away some super-badly-tuned cells just to denoise a bits
%   note: currently, the lower 95% C.I. of R2 must be above this value
pars.min_r2 = 0;
pars.n_bins_downs = 3;

% simple tuning script
pars.bad_neuron_params.min_fr           = 0.1;
% Do shunt check
pars.bad_neuron_params.shunt_check_yn   = false;

% I did 1000 for the saved figures, but 100 is good enough in general
pars.tuning_conf_params = {'bootstrap',100,0.95};

switch lower(pars.monkey)
    case {'chips'} % S1
        pars.mani_dims = 1:8; % Neural modes to use
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 19};
    case {'han'} % S1
        pars.mani_dims = 1:8;  % Neural modes to use
        pars.idx_start          = {'idx_go_cue', 0};
        pars.idx_end            = {'idx_go_cue', 22};
    case {'chewie','chewie2','mihili','mrt','jaco'}
        switch pars.spiking_inputs{1}
            case 'M1_spikes'
                pars.mani_dims = 1:10; % Neural modes to use
%                 pars.idx_start  = {'idx_movement_on',-2}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
%                 pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
                pars.idx_start  = {'idx_movement_on',-4}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
                pars.idx_end    = {'idx_movement_on',13}; % {'idx_movement_on',13}% {''}; % {'idx_go_cue',18}
            case 'PMd_spikes'
                pars.mani_dims = 1:16; % Neural modes to use
                pars.idx_start  = {'idx_go_cue',-13};
                pars.idx_end    = {'idx_go_cue',2};
        end
end

monkey_file_lists;


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
    elseif strcmpi(pars.monkey,'han')
        date = datestr(datenum(temp{3},'yyyy-mm-dd'),'mm-dd-yyyy');
        task = 'COactpas';
        monkey = 'Han';
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
if strcmpi(pars.monkey,'han') || strcmpi(pars.monkey,'chips')
    load_td = [];
    for iFile = 1:length(file_list)
        load(fullfile(data_dir,pars.monkey,file_list{iFile}));
        [~,td] = getTDidx(trial_data,'result','R');
        td = stripSpikeSorting(td);
        td = removeBadTrials(td,struct('nan_idx_names','idx_go_cue'));
        td = binTD(td,pars.n_bins_downs);
        td = removeBadNeurons(td,pars.bad_neuron_params);
        td = trimTD(td,pars.idx_start,pars.idx_end);
        
        load_td = [load_td, td];
    end
    trial_data = load_td; clear load_td td;
    
else % chewie and mihili both fall in this one
    filepaths = cell(size(file_list));
    for iFile = 1:length(file_list)
        filepaths{iFile} = fullfile(data_dir,pars.monkey,file_list{iFile});
    end
    load_fn_call = { ...
        {@getTDidx,'result','R','epoch','BL'}, ...
        @stripSpikeSorting, ...
        @removeBadTrials, ...
        {@binTD,pars.n_bins_downs}, ...
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

% figure('Position',[100 100 200 600]);
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
    
%     if iFile == 8
%         idx = find(td(1).M1_unit_guide(:,1)==28);
%         th = -pi:pi/32:pi;
%         
%         subplot(2,1,1); hold all;
%         plot(dir*180/pi,fr(:,idx),'.','markersize',16,'color',[0.3 0.3 0.3]);
%         plot(th*180/pi,tc(idx,1)+tc(idx,2)*cos(th - tc(idx,3)),'r-','LineWidth',2);
%         set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-180 180],'YLim',[0 300]);
%         set(gca,'XTick',[-90 0 90],'XTickLabel',[-90 0 90]);
%         xlabel('Reach Direction');
%         ylabel('Firing Rate');
%     elseif iFile == 14
%         idx = find(td(1).M1_unit_guide(:,1)==28);
%         th = -pi:pi/32:pi;
%         
%         subplot(2,1,2); hold all;
%         plot(dir*180/pi,fr(:,idx),'.','markersize',16,'color',[0.3 0.3 0.3]);
%         plot(th*180/pi,tc(idx,1)+tc(idx,2)*cos(th - tc(idx,3)),'r-','LineWidth',2);
%         set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-180 180],'YLim',[0 300]);
%         set(gca,'XTick',[-90 0 90],'XTickLabel',[-90 0 90]);
%         xlabel('Reach Direction');
%         ylabel('Firing Rate');
%     end
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
    
    fn = fullfile(save_dir,'Tuning stability',['NeuralTuningStability_' pars.monkey '_' pars.array '.mat']);
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
        
        plot(x,m,'.','markersize',32,'color',[0.3 0.3 0.3]);
        %plot([x, x],[m-s,m+s],'-','LineWidth',2,'color',[0.3 0.3 0.3])
        
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
    fig1=gcf;
fig1.Renderer='Painters';
    fn = fullfile(save_dir,'Tuning stability',['NeuralTuningStability_' pars.monkey '_' pars.array]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end




%% now plot stuff

figure('Position',[100 100 800 400]);
for tp = 1:length(which_tuning_params)
    
    subplot(1,length(which_tuning_params),tp);
    hold all;
    
    [all_diff,all_days] = deal([]);
    for c = 1:n_comb_sessions
        m = tuneDiff_m(c,tp);
%         s = tuneDiff_sd(c,tp);
        d  = tuneDiff_all{c,tp};
        all_diff = [all_diff; d];
        
        x = datenum(file_info(comb_sessions(c,2)).date,'mm-dd-yyyy') - datenum(file_info(comb_sessions(c,1)).date,'mm-dd-yyyy');
        all_days = [all_days; repmat(x,size(d))];
        
        plot(x,m,'.','markersize',32,'color',[0.3 0.3 0.3]);
%         plot([x, x],[m-s,m+s],'-','LineWidth',2,'color',[0.3 0.3 0.3])
        
    end
    
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(diff_days)+1]);
    
    % fit a line
    [b,bint,~,~,s] = regress(all_diff,[ones(size(all_diff)) all_days]);
    r = s(1);
    p = s(3);
    
    V = axis;
    h = patch([ ...
        V(1):V(2) ...
        fliplr(V(1):V(2))], ...
        [bint(1,1) + bint(2,1)*(V(1):V(2)) ...
        fliplr(bint(1,2) + bint(2,2)*(V(1):V(2)))], ...
        'k--','LineWidth',2);
    set(h,'EdgeColor','none','EdgeAlpha',0,'FaceAlpha',0.2);
    plot(V(1:2),b(1)+b(2)*V(1:2),'k-','LineWidth',2);

    text(0.1*V(2),0.1*V(3),['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x; R^2 = ' num2str(r,3) '; p = ' num2str(p,3)],'FontSize',12);
    
    title([pars.monkey ' - ' pars.array]);
    xlabel('Days between sessions')
    ylabel(['|Change in ' pars.tuning_param_labels{tp} '|']);
    
end

if save_figs
        fig1=gcf;
fig1.Renderer='Painters';
    fn = fullfile(save_dir,'Tuning stability',['NeuralTuningStability_' pars.monkey '_' pars.array '_fiterr']);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end



