%
%
% Compare poulation dynamics across days using the TD struct
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params and files

clear all, close all

% day1        = '09-09-2016';
arrays      = {'M1_array'}; % This is yet to be used

idx_start   = {'idx_go_cue',0};
idx_end     = {'idx_go_cue',60};

mani_dims   = 1:10;

% files       = { ...
%                 'Chewie_CO_VR_2013-10-03.mat', ...
%                 'Chewie_CO_FF_2013-10-22.mat', ...
%                 'Chewie_CO_FF_2013-10-23.mat', ...
%                 'Chewie_CO_FF_2013-10-31.mat', ...
%                 'Chewie_CO_FF_2013-11-01.mat', ...
%                 'Chewie_CO_FF_2013-12-03.mat', ...
%                 'Chewie_CO_FF_2013-12-04.mat', ...
%                 'Chewie_CO_VR_2013-12-19.mat', ...
%                 'Chewie_CO_VR_2013-12-20.mat'
%                 };

% files       = { ...
%                 'Chewie_CO_FF_2015-06-29.mat', ...
%                 'Chewie_CO_FF_2015-06-30.mat', ...
%                 'Chewie_CO_FF_2015-07-01.mat', ...
%                 'Chewie_CO_FF_2015-07-03.mat', ...
%                 'Chewie_CO_FF_2015-07-06.mat', ...
%                 'Chewie_CO_FF_2015-07-07.mat', ...
%                 'Chewie_CO_FF_2015-07-08.mat', ...
%                 'Chewie_CO_VR_2015-07-09.mat', ...
%                 'Chewie_CO_VR_2015-07-10.mat', ...
%                 'Chewie_CO_VR_2015-07-13.mat', ...
%                 'Chewie_CO_VR_2015-07-14.mat', ...
%                 'Chewie_CO_VR_2015-07-15.mat', ...
%                 'Chewie_CO_VR_2015-07-15.mat' ...
%                 };
                

files       = { ...
                'Chewie_CO_VR_2016-09-09.mat', ...
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

% files           = { ...
%                     'Mihili_CO_FF_2014-02-03.mat', ...
%                     'Mihili_CO_FF_2014-02-17.mat', ...
%                     'Mihili_CO_FF_2014-02-18.mat', ...
%                     'Mihili_CO_VR_2014-03-03.mat', ...
%                     'Mihili_CO_VR_2014-03-04.mat', ...
%                     'Mihili_CO_VR_2014-03-06.mat', ...
%                     'Mihili_CO_FF_2014-03-07.mat' ...
%                     };
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data, unsorting the channels. Baseline only

master_td   = loadTDfiles(  files, ...
                            {@mergeSortedNeurons}, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@trimTD,idx_start,idx_end}, ...
                            {@sqrtTransform,{'M1_spikes'}}, ...
                            {@smoothSignals,struct('signals',{'M1_spikes'},'calc_fr',true,'kernel_SD',0.05)}, ...
                            {@getPCA,struct('signals',{'M1_spikes'})} ...
                         );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get some meta information

% get the days
sessions    = unique({master_td.date});
% the targets
targets     = unique(cell2mat({master_td.target_direction}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the channels that we have in each day


% elecs       = zeros(length(sessions),96);
% 
% for s = 1:length(sessions)
%     % find the index of the first trial of each session
%     idx1    = cell2mat( arrayfun( @(x) strncmp( x.date, sessions{s}, length(sessions{s}) ), ...
%                 master_td, 'UniformOutput', false ) );
%     idx1    = find(idx1,1);
%     % store the electrodes in this session
%     elecs(s,unique(master_td(idx1).M1_unit_guide)) = 1;
% end
% 
% figure,imagesc(~elecs),colormap('bone'),
% set(gca,'TickDir','out','FontSize',14), box off
% xlabel('Electrode'),ylabel('Session'),title('Electrodes with neurons')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalize number of trials across sessions and targets

master_td   = equalNbrTrialsSessions( master_td );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do CCA across sessions

% get all pairwise comparisons of sessions
comb_sessions   = nchoosek(1:length(sessions),2);

% do
for c = 1:size(comb_sessions,1)
    
    % get the data
    trials1 = getTDidx( master_td, 'date', sessions{comb_sessions(c,1)},...
                    'target_direction', targets(1:length(targets)) );
    trials2 = getTDidx( master_td, 'date', sessions{comb_sessions(c,2)},...
                    'target_direction', targets(1:length(targets)) );
    
    %trials2 = trials2(randperm(length(trials2)));
                
    % compare dynamics with CCA
    cca_info(c) = compDynamics( master_td, 'M1_pca', trials1, trials2, mani_dims );
    
    % compare dynamics with good old forrelations
    corr_info(c) = corrDynamics( master_td, 'M1_pca', trials1, trials2, mani_dims );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do CCA within sessions to have a "baseline condition"

for s = 1:length(sessions)
   
    % retrieve session
    [trialss, tds] = getTDidx( master_td, 'date', sessions{s},...
                    'target_direction', targets(1:length(targets)) );
    
    trialsptarget = (length(trials1)/length(targets));
    htrials = floor(trialsptarget/2);
    h1      = 1:htrials;
    h2      = (trialsptarget-htrials+1):trialsptarget;
                
    % get first and second half of the trials
    trials1h = [];
    trials2h = [];
    for t = 1:length(targets)
        trials1h = [trials1h h1+(t-1)*trialsptarget];
        trials2h = [trials2h h2+(t-1)*trialsptarget];
    end
               
    % compare dynamics with CCA
    cca_within(s) = compDynamics( tds, 'M1_pca', trials1h, trials2h, mani_dims );
    
    % compare dynamics with good old forrelations
    corr_within(s) = corrDynamics( tds, 'M1_pca', trials1h, trials2h, mani_dims );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get time between sessions

diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pool all CCs
all_ccs     = cell2mat({cca_info.cc}');
all_corrs   = cell2mat({corr_info.r}');

% summary statistics
dims_stats  = 1:3;

mn_cc           = zeros(size(comb_sessions,1),1);
sem_cc          = zeros(size(comb_sessions,1),1);
mn_corr         = zeros(size(comb_sessions,1),1);
sem_corr        = zeros(size(comb_sessions,1),1);


for i = 1:size(comb_sessions,1)
    cc          = abs(all_ccs(i,dims_stats));
    cr          = abs(all_corrs(i,dims_stats));
    mn_cc(i)    = mean(cc);
    sem_cc(i)   = std(cc)/sqrt(numel(mani_dims));
    mn_corr(i)  = mean(cr);
    sem_corr(i) = std(cr)/sqrt(numel(mani_dims));
end

% linear fits
fit_cc          = polyfit(diff_days,mn_cc',1);
fit_corr        = polyfit(diff_days,mn_corr',1);
x_plots         = [min(diff_days)-1, max(diff_days)+1];
y_cc            = polyval(fit_cc,x_plots);
y_corr          = polyval(fit_corr,x_plots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairwise corrs & Canon corrs vs dimensions

figure,plot(all_ccs'),hold on, plot(abs(all_corrs),'k')
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Projection'),ylabel('Correlation')
ylim([0 1]),xlim([min(mani_dims), max(mani_dims)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairwise corrs & Canon corrs vs dimensions
figure, hold on
errorbar(diff_days,mn_cc,sem_cc,'.k','markersize',32)
errorbar(diff_days,mn_corr,sem_corr,'o','color',[.6 .6 .6],'markersize',10)
plot(x_plots,y_cc,'k','linewidth',1.5)
plot(x_plots,y_corr,'color',[.6 .6 .6],'linewidth',1.5)
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days btw sessions'),ylabel('Correlation')
ylim([0 1]), xlim([0 max(diff_days)+1])
legend('aligned','unaligned','Location','SouthEast'),legend boxoff


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairwise corrs & Canon corrs vs dimensions
figure, hold on
errorbar(diff_days,mn_cc,sem_cc,'.k','markersize',32)
errorbar(diff_days,mn_corr,sem_corr,'o','color',[.6 .6 .6],'markersize',10)
plot(x_plots,y_cc,'k','linewidth',1.5)
plot(x_plots,y_corr,'color',[.6 .6 .6],'linewidth',1.5)
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days btw sessions'),ylabel('Correlation')
ylim([0 1]), xlim([0 max(diff_days)+1])
legend('aligned','unaligned','Location','SouthEast'),legend boxoff

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pool all within day CCs
all_within_ccs = cell2mat({cca_within.cc}');
all_within_corrs = cell2mat({corr_within.r}');

% summary statistics --taking the same number of dimenions as above

mn_wcc          = zeros(length(sessions),1);
sem_wcc         = zeros(length(sessions),1);
mn_wcorr        = zeros(length(sessions),1);
sem_wcorr       = zeros(length(sessions),1);


for i = 1:length(sessions)
    wcc         = abs(all_within_ccs(i,dims_stats));
    wcr         = abs(all_within_corrs(i,dims_stats));
    mn_wcc(i)   = mean(wcc);
    sem_wcc(i)  = std(wcc)/sqrt(numel(mani_dims));
    mn_wcorr(i) = mean(wcr);
    sem_wcorr(i) = std(wcr)/sqrt(numel(mani_dims));
end

% append number of within sessions
diff_days_all   = [diff_days zeros(1,length(sessions))];

% append correlations
mn_cc_all       = [mn_cc; mn_wcc]';
sem_cc_all      = [sem_cc; sem_wcc]';
mn_corr_all     = [mn_corr; mn_wcorr]';
sem_corr_all    = [sem_corr; sem_wcorr]';


% linear fits including these baseline days
fit_cc_all      = polyfit(diff_days_all,mn_cc_all,1);
fit_corr_all    = polyfit(diff_days_all,mn_corr_all,1);
x_plots_all     = [min(diff_days_all)-2, max(diff_days_all)+1];
y_cc_all        = polyval(fit_cc_all,x_plots_all);
y_corr_all      = polyval(fit_corr_all,x_plots_all);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairwise corrs & Canon corrs vs dimensions including basline
figure, hold on
errorbar(diff_days_all,mn_cc_all,sem_cc_all,'.k','markersize',32)
errorbar(diff_days_all,mn_corr_all,sem_corr_all,'o','color',[.6 .6 .6],'markersize',10)
plot(x_plots_all,y_cc_all,'k','linewidth',1.5)
plot(x_plots_all,y_corr_all,'color',[.6 .6 .6],'linewidth',1.5)
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days btw sessions'),ylabel('Correlation')
ylim([0 1]), xlim([-1 max(diff_days_all)+1])
legend('aligned','unaligned','Location','SouthEast'),legend boxoff
