clear;
close all;
clc;

[save_dir, data_dir] = get_computer_paths();
save_them = false;

pars.array          = 'PMd';
pars.spiking_inputs = {[pars.array '_spikes']};


switch pars.array
    case {'M1','PMd'}
        pars.monkey         = 'Chewie'; % 'chewie2'; 'chewie'; 'mihili'; 'han'; 'chips'; 'jaco'
        
        fr_max = 180;
vel_lim = [-25 25];

        date0 = '2016-09-09';
        date1 = '2016-10-05';
        date2 = '2016-10-21';
        
        tdfile1 = ['/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_CO_FF_' date1 '.mat'];
        tdfile2 = ['/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_CO_CS_' date2 '.mat'];
        
        %--------------------------------------------------------------------------
        % load CDS for first day
        load(['/Users/mattperich/Desktop/Chewie/CDS/' date1 '/Chewie_CO_FF_BL_10052016_001.mat']);
        cds1 = cds;
        % load CDS for second day
        load(['/Users/mattperich/Desktop/Chewie/CDS/' date2 '/Chewie_CO_CS_BL_10212016_001.mat']);
        cds2 = cds;
        clear cds;
        
        params_stability_over_time;
        
        [master_td, pars_td] = loadTDfiles(  {tdfile1,tdfile2}, ...
            {@getTDidx,'epoch','BL','result','R'}, ...    {@stripSpikeSorting},...
            {@binTD,pars.n_bins_downs}, ...
            {@removeBadNeurons,pars.bad_neuron_params},...
            {@sqrtTransform,pars.spiking_inputs}, ...
            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
            {@getPCA,struct('signals',pars.spiking_inputs)}, ...
            {@trimTD,pars.idx_start,pars.idx_end});
        
        master_td = removeBadTrials(master_td,struct('nan_idx_names','idx_movement_on'));
        
    case 'S1'
        pars.monkey = 'Chips';
        
        fr_max = 180;
vel_lim = [-30 30];

        date0 = '2015-11-13';
        date1 = '2015-11-13';
        date2 = '2015-12-11';
        
        tdfile1 = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chips/Chips_20151113_TD_nosort_notrack_noemg.mat';
        tdfile2 = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chips/Chips_20151211_TD_nosort_notrack_noemg.mat';
        
        params_stability_over_time;
        
        [master_td, pars_td] = loadTDfiles(  {tdfile1,tdfile2}, ...
            {@getTDidx,'result','R'}, ...    {@stripSpikeSorting},...
            {@binTD,pars.n_bins_downs}, ...
            {@removeBadNeurons,pars.bad_neuron_params},...
            {@removeBadTrials,struct('nan_idx_names','idx_go_cue','ranges',{{'idx_go_cue','idx_trial_end',[19 Inf]}})}, ...
            {@sqrtTransform,pars.spiking_inputs}, ...
            {@smoothSignals,struct('signals',pars.spiking_inputs,'calc_fr',true,'kernel_SD',pars.kernel_SD)}, ...
            {@trimTD,{'idx_target_on',0},{'idx_trial_end',0}}, ...
            {@getPCA,struct('signals',pars.spiking_inputs)}, ...
            {@trimTD,pars});
end

%--------------------------------------------------------------------------
% load TD for both days

% some sessions are missing force so just throw it out
% since we don't n eed it for this analysis
master_td = rmfield(master_td,'force');
[master_td, num_trials] = equalNbrTrialsSessions(master_td);
master_td = trimTD(master_td,pars);

% get the trial data and unit idx belonging to this electrode
[~,td1] = getTDidx(master_td,'date',datestr(datenum(date1,'yyyy-mm-dd'),'mm-dd-yyyy'));
[~,td2] = getTDidx(master_td,'date',datestr(datenum(date2,'yyyy-mm-dd'),'mm-dd-yyyy'));

d_day1 = datenum(date1) - datenum(date0);
d_day2 = datenum(date2) - datenum(date0);


%% make a proper raster plot for the two days
do_norm = false;
n_whitespace = 1;

%--------------------------------------------------------------------------
close all;

td1_temp = td1;
td2_temp = td2;

% td1_temp = stripSpikeSorting(td1_temp);
% td2_temp = stripSpikeSorting(td2_temp);

td1_temp = trialAverage(td1_temp,'target_direction');
td2_temp = trialAverage(td2_temp,'target_direction');

% add some whitespace to each trial
for trial = 1:length(td1_temp)
    td1_temp(trial).([pars.array '_spikes'])  = cat(1,td1_temp(trial).([pars.array '_spikes']),NaN(n_whitespace,size(td1_temp(trial).([pars.array '_spikes']),2)));
    td1_temp(trial).vel  = cat(1,td1_temp(trial).vel,NaN(n_whitespace,size(td1_temp(trial).vel,2)));
end
% add some whitespace to each trial
for trial = 1:length(td2_temp)
    td2_temp(trial).([pars.array '_spikes'])  = cat(1,td2_temp(trial).([pars.array '_spikes']),NaN(n_whitespace,size(td2_temp(trial).([pars.array '_spikes']),2)));
    td2_temp(trial).vel  = cat(1,td2_temp(trial).vel,NaN(n_whitespace,size(td2_temp(trial).vel,2)));
end

figure('Position',[100 100 1000 500]);

%--------------------------------------------------------------------------
[~,idx]  = get_test_train_trials(td1_temp,1);
fr = cat(1,td1_temp(idx).([pars.array '_spikes']))./0.01;
if do_norm
    fr = fr./repmat(nanmean(fr,1),size(fr,1),1);
end
fr(isnan(fr)) = Inf;
size(fr)
% plot spikes
ax(1) = subplot(4,2,[1,3,5]);
imagesc(fr');
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'CLim',[0 fr_max]);
axis tight;
ylabel('Sorted neurons')
title(['Day ' num2str(d_day1)]);
colormap hot;
% colormap(brewermap(100,'*Blues'));

% now plot vel
ax(3) = subplot(4,2,7);
plot(cat(1,td1_temp(idx).vel),'LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YLim',vel_lim);
xlabel('Time (10ms bins)');

%--------------------------------------------------------------------------
[~,idx]  = get_test_train_trials(td2_temp,1);
fr = cat(1,td2_temp(idx).([pars.array '_spikes']))./0.01;
if do_norm
    fr = fr./repmat(nanmean(fr,1),size(fr,1),1);
end
fr(isnan(fr)) = Inf;

% plot spikes
ax(2) = subplot(4,2,[2,4,6]);
imagesc(fr');
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'CLim',[0 fr_max]);
axis tight;
title(['Day ' num2str(d_day2)]);

% now plot vel
ax(4) = subplot(4,2,8);
plot(cat(1,td2_temp(idx).vel),'LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YLim',vel_lim);

if save_them
    fn = fullfile(save_dir,'Neural activity',[pars.monkey '_' pars.array '_PopulationRasters']);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end

%% pick representative electrode(s) that have neurons on both days
elec = 32; % 13, 28, 32, 50, 92

% some parameters
fr_max = 30;
wf_lims = [-600, 500];
isi_bins = 0:250;

%--------------------------------------------------------------------------
% get the waveforms for the electrode
idx1 = find(td1(1).M1_unit_guide(:,1) == elec);
idx2 = find(td2(1).M1_unit_guide(:,1) == elec);

wf1 = cell(1,length(idx1));
for i = 1:length(idx1)
    idx = [cds1.units.chan] == td1(1).M1_unit_guide(idx1(i),1) & ...
        [cds1.units.ID] == td1(1).M1_unit_guide(idx1(i),2) & ...
        strcmpi({cds1.units.array},'M1');
    wf1{i} = cat(1,cds1.units(idx).spikes.wave);
end

wf2 = cell(1,length(idx2));
for i = 1:length(idx2)
    idx = [cds2.units.chan] == td2(1).M1_unit_guide(idx2(i),1) & ...
        [cds2.units.ID] == td2(1).M1_unit_guide(idx2(i),2) & ...
        strcmpi({cds2.units.array},'M1');
    wf2{i} = cat(1,cds2.units(idx).spikes.wave);
end

%--------------------------------------------------------------------------
% get the ISIs for the electrode
isi1 = zeros(length(idx1),length(isi_bins)-1);
for i = 1:length(idx1)
    idx = [cds1.units.chan] == td1(1).M1_unit_guide(idx1(i),1) & ...
        [cds1.units.ID] == td1(1).M1_unit_guide(idx1(i),2) & ...
        strcmpi({cds1.units.array},'M1');
    isi = diff(cds1.units(idx).spikes.ts);
    isi(isi > 0.25) = [];
    
    isi1(i,:) = histcounts(isi*1000,isi_bins);
end

isi2 = zeros(length(idx2),length(isi_bins)-1);
for i = 1:length(idx2)
    idx = [cds2.units.chan] == td2(1).M1_unit_guide(idx2(i),1) & ...
        [cds2.units.ID] == td2(1).M1_unit_guide(idx2(i),2) & ...
        strcmpi({cds2.units.array},'M1');
    isi = diff(cds2.units(idx).spikes.ts);
    isi(isi > 0.25) = [];
    
    isi2(i,:) = histcounts(isi*1000,isi_bins);
end

% normalize a bit for plotting
isi1 = isi1./repmat(sum(isi1,2),1,size(isi1,2));
isi2 = isi2./repmat(sum(isi2,2),1,size(isi2,2));

%--------------------------------------------------------------------------
% now do the plotting
figure('Position',[100 100 800 800]);
subplot1(5,5);
% plot the waveforms
subplot1(1); hold all;
for i = 1:length(idx1)
    plot(mean(wf1{i},1),'LineWidth',2)
end
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[1 48],'YLim',wf_lims);
title(['Day ' num2str(d_day1)]);

%--------------------------------------------------------------------------
% plot the ISIs
subplot1(5); hold all;
plot(isi1','LineWidth',2)
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 isi_bins(end)],'YLim',[0 0.05]);
set(gca,'XTick',isi_bins(1:100:end),'XTickLabel',isi_bins(1:100:end));
xlabel('ISI');

%--------------------------------------------------------------------------
% plot the PETHs
plot_idx = [17,23,19,15,9,3,7,11];
utheta = unique([td1.target_direction]);

fr_dir = [];
for u = 1:length(utheta)
    [~,td] = getTDidx(td1,'target_direction',utheta(u));
    
    temp = cat(3,td.M1_spikes);
    temp = temp(:,idx1,:);
    temp = squeeze(sum(temp,3));
    
    fr_dir = cat(3,fr_dir,temp);
    
    subplot1(plot_idx(u));
    hold all;
    plot(squeeze(fr_dir(:,:,u)),'LineWidth',2);
    plot(repmat(abs(pars.idx_start{2}),1,2),[0 fr_max],'k--')
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0 fr_max]);
    axis square;
    
    if u ~= 2
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    else
        xlabel('10 ms bins');
        ylabel(['# spikes (' num2str(num_trials) ' trials)']);
        set(gca,'YTick',0:10:fr_max,'YTickLabel',0:10:fr_max);
    end
    
    if u == 6
        title(['Electrode ' num2str(elec)]);
    end
end
for i = 1:25
    subplot1(i);
    if ~ismember(i,plot_idx) && i ~= 1 && i ~= 5
        set(gca,'Visible','off');
    end
end


if save_them
    fn = fullfile(save_dir,'Neural activity',[pars.monkey '_' pars.array '_elec' num2str(elec) '_day' num2str(d_day1)]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end




% now plot the second day
figure('Position',[100 100 800 800]);
subplot1(5,5);
%--------------------------------------------------------------------------
% plot the waveforms
subplot1(1); hold all;
for i = 1:length(idx2)
    plot(mean(wf2{i},1),'LineWidth',2)
end
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[1 48],'YLim',wf_lims);
title([ 'Day ' num2str(d_day2) ]);

%--------------------------------------------------------------------------
% plot the ISIs
subplot1(5); hold all;
plot(isi2','LineWidth',2)
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 isi_bins(end)],'YLim',[0 0.05]);
set(gca,'XTick',isi_bins(1:100:end),'XTickLabel',isi_bins(1:100:end));
xlabel('ISI (s)');


%--------------------------------------------------------------------------
fr_dir = [];
for u = 1:length(utheta)
    [~,td] = getTDidx(td2,'target_direction',utheta(u));
    
    temp = cat(3,td.M1_spikes);
    temp = temp(:,idx2,:);
    temp = squeeze(sum(temp,3));
    
    fr_dir = cat(3,fr_dir,temp);
    
    subplot1(plot_idx(u));
    hold all;
    plot(squeeze(fr_dir(:,:,u)),'LineWidth',2);
    plot(repmat(abs(pars.idx_start{2}),1,2),[0 fr_max],'k--');
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0 fr_max]);
    axis square;
    
    if u ~= 2
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    else
        xlabel('10 ms bins');
        ylabel(['# spikes (' num2str(num_trials) ' trials)']);
        set(gca,'YTick',0:10:fr_max,'YTickLabel',0:10:fr_max);
    end
    
    if u == 6
        title(['Electrode ' num2str(elec)]);
    end
end

for i = 1:25
    subplot1(i);
    if ~ismember(i,plot_idx) && i ~= 1 && i ~= 5
        set(gca,'Visible','off');
    end
end

if save_them
    fn = fullfile(save_dir,'Neural activity',[pars.monkey '_' pars.array '_elec' num2str(elec) '_day' num2str(d_day2)]);
    savefig(gcf,[fn '.fig']);
    saveas(gcf,fn,'png');
    saveas(gcf,fn,'pdf');
end
