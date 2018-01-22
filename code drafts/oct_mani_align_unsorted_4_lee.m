
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

% Parameters / definition
monkey          = 'chewie'; % 'mihili'; 'chewie'

unsorted_yn     = true;

spiking_inputs  = {'M1_spikes'}; % 'PMd_spikes'

% The behavioral variable to be predicted
dec_var         = 'vel';
% Decoder inputs
dec_input       = 'aligned_data'; % 'aligned_data'; 'raw_data'

% When we start looking at each trial
idx_start       = {'idx_go_cue',0}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
% When we stop looking at each trial --if empty, the duration of the shortest trial 
idx_end         = {'idx_go_cue',50}; % {''}; % {'idx_go_cue',10}
% Gaussian kernel for smooting
kernel_SD       = 0.1;
% Neural modes that will be aligned
mani_dims       = 1:10;

% Downsampling rate (nbr of bins that will be combined)
n_bins          = 5;
% N bins for decoder history
hist_bins       = 3;
% N folds for decoder cross-validation
n_folds         = 6; % this gives blocks of 20 trials for chewie

% Minimum FR per channel
min_fr          = 0.1;

% Flag to only keep electrodes that have units in all sessions
only_electrodes_w_units_all_sessions = false;

% Choose sesssions to discard if any
sessions_discard = [];


master_td       = loadTDfiles(  '/Users/juangallego/Documents/NeuroPlast/Data/Chewie/Chewie_unsorted_M1_only_stab_over_time_Ali.mat', ...
                            {@getTDidx,'result','R'});

                        
                        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get some meta info about the data

% get the sessions
sessions    = unique({master_td.date});

% get the targets
targets     = unique(cell2mat({master_td.target_direction}));




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only keep a subset of sessions, to avoid having too few sessions


discard = getTDidx(master_td,{'date',{sessions{sessions_discard}}});
master_td_orig = master_td;

% discard sessions 
master_td(discard) = [];
sessions(sessions_discard) = [];




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim the trials

% get the minimum trial duration, if needed for idx_end
if isempty(idx_end{1})
    min_mov_duration = min( arrayfun( @(x) x.idx_trial_end-x.idx_movement_on, master_td) );
    idx_end = {'idx_movement_on',min_mov_duration};
end    

master_td   = trimTD( master_td, idx_start, idx_end );




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some pre-processing

% bin the data (downsample)
for s = 1:length(sessions)
    this_s = getTDidx(master_td,{'date',sessions{s}});
    master_td(this_s) = binTD(master_td(this_s),n_bins);
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep good channels only: delete channels that are not units, and get rid
% of shunted electrodes, or very low firing channels (optional)


n_units     = zeros(1,length(sessions));
for s = 1:length(sessions)
    this_s = getTDidx(master_td,{'date',sessions{s}});
    
    % ---------------------------------------------------------------------
    % 1. Delete channels >96 (the darned ch 137)
    [~,idx_no_uea] = setdiff(master_td(this_s(1)).M1_unit_guide(:,1),1:96);
    

    % ---------------------------------------------------------------------
    % 2. Find shunted channels for each session
    [~,sts] = removeBadNeurons(master_td(this_s),struct('do_shunt_check',true));
    shunted_session{s} = sts;

    
    % ---------------------------------------------------------------------
    % 3. Find channels with 0 spikes
    [~,zspks] = removeBadNeurons(master_td(this_s),struct('do_fr_check',true,'min_fr',min_fr));
    zero_spike_chs_session{s} = zspks;

    
    % ---------------------------------------------------------------------
    % 4. Delete
    
    idx_del = unique([idx_no_uea, sts, zspks]);
    
    for t = 1:length(this_s)
        master_td(this_s(t)).M1_unit_guide(idx_del,:)=[];
        master_td(this_s(t)).M1_spikes(:,idx_del)=[];
    end
    
    % ---------------------------------------------------------------------
    % store the nbr of good units per session
    n_units(s) = size(master_td(this_s(1)).M1_unit_guide,1);
    % store channel numbers good units
    units_session{s} = master_td(this_s(1)).M1_unit_guide(:,1);
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get electrodes that are usable across all sessions

% NOTE: THIS may be the same as Matt's function getCommonUnits !!!!

% array_units_session = zeros(length(sessions),96);
% for s = 1:length(sessions)
%     array_units_session(s,:) = ismember(1:96,units_session{s});
% end
% saus = sum(array_units_session,1);
% usable_all = saus==size(array_units_session,1);
% 
% % plot
% figure,subplot(211)
% imagesc(~array_units_session),colormap('gray')
% set(gca,'TickDir','out','FontSize',14), ylabel('Session number')
% title(['Min FR = ' num2str(min_fr) ' - Nbr. common electrodes = ' num2str(sum(usable_all))])
% subplot(212)
% imagesc(~usable_all)
% set(gca,'TickDir','out','FontSize',14), xlabel('Electrode number')
% set(gca,'YTick',[]),ylabel('Electrodes all sessions')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of electrodes that are not usable across all sessions

if only_electrodes_w_units_all_sessions
    
    uea_chs = 1:96;
    usable_uea_chs = uea_chs(usable_all);

    for s = 1:length(sessions)

        this_s = getTDidx(master_td,{'date',sessions{s}});

        % Find index electrodes to delete
        this_ch = master_td(this_s(1)).M1_unit_guide(:,1);
        [~,idx_del] = setdiff(this_ch,usable_uea_chs);

        % Delete these electrodes
        for t = 1:length(this_s)
            master_td(this_s(t)).M1_unit_guide(idx_del,:)=[];
            master_td(this_s(t)).M1_spikes(:,idx_del)=[];
        end
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep the same number of trials for all the targets and sessions --this is
% important to align the neural mode dynamics with CCA

master_td   = equalNbrTrialsSessions(master_td);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do PCA --only necessary for unsorted datafiles, because in this case we
% don't do PCA during the loading (because we want to have the same
% electrodes in all the sessions)

% if unsorted_yn
%     
%     for s = 1:length(sessions)
%         this_s = getTDidx( master_td,{'date',sessions{s}} );
%         [master_td_new(this_s), pca_info(s)] = getPCA( master_td(this_s), struct(...
%                                     'signals',spiking_inputs,...
%                                     'sqrt_transform',true,...
%                                     'do_smoothing',true,...
%                                     'kernelSD',0.1 ));
%     end
%     
%     master_td = master_td_new;
% end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do CCA of the neural trajectories across sessions, to align the dynamics
% and see how similar they are

% % get all pairwise comparisons of sessions
% comb_sessions   = nchoosek(1:length(sessions),2);
% 
% % do
% disp('Aligning neural mode dynamics across days')
% for c = 1:size(comb_sessions,1)
%     
%     % get the data
%     trials1 = getTDidx( master_td, 'date', sessions{comb_sessions(c,1)} );
%     trials2 = getTDidx( master_td, 'date', sessions{comb_sessions(c,2)} );
%     
%     %trials2 = trials2(randperm(length(trials2)));
% 
%     switch spiking_inputs{1}
%         case 'M1_spikes'
%             manifold = 'M1_pca';
%         case 'PMd_spikes'
%             manifold = 'PMd_pca';
%     end
%     % compare dynamics with CCA
%     cca_info(c) = compDynamics( master_td, manifold, trials1, trials2, mani_dims );
%     
%     % compare dynamics with good old forrelations
%     corr_info(c) = corrDynamics( master_td, manifold, trials1, trials2, mani_dims );
% end
% 
% 
% % get number of days between sessions, to summarize the results
% diff_days   = zeros(1,size(comb_sessions,1));
% for c = 1:size(comb_sessions,1)
%     diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
% end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build a decoder on one day using threshold crossings and test it in another

% get all pairwise comparisons of sessions
comb_sessions   = nchoosek(1:length(sessions),2);

% do
disp('Building decoders and testing them in another session')
for c = 1:size(comb_sessions,1)
    
    % get the data
    trials1 = getTDidx( master_td, 'date', sessions{comb_sessions(c,1)} );
    trials2 = getTDidx( master_td, 'date', sessions{comb_sessions(c,2)} );

    % get the channels that have threshold crossings in both datasets
    this_td = master_td([trials1 trials2]);
    this_td = getCommonUnits(this_td);
    
    % split the build and test data, and duplicate and shift
    [~,td1] = getTDidx(this_td, 'date', sessions{comb_sessions(c,1)} );
    [~,td2] = getTDidx(this_td, 'date', sessions{comb_sessions(c,2)} );
    
    td1 = dupeAndShift(td1,{'M1_spikes',hist_bins});
    td2 = dupeAndShift(td2,{'M1_spikes',hist_bins});
    
    % define model parameters --move up!!!
    mod_params.in_signals = 'M1_spikes_shift';
    mod_params.model_type = 'linmodel';
    mod_params.out_signals = dec_var;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build the model on dataset 1

    [td1, mod_info] = getModel(td1,mod_params);

%     % Test model fit without cross-validation
%     % compute the R2 of the within day predictions for the first date
%     Y1 = get_vars(td1,{dec_var,1:size(td1(1).(dec_var),2)});
%     Y1hat = get_vars(td1,{'linmodel_default',1:size(td1(1).(dec_var),2)});
% 
%     R2 = compute_r2(Y1,Y1hat,'corr');
%     withinR2(c,:) = R2;

    % Do multi-fold cross-validation
    % -- randomize trial order, because now they are sorted by target, and
    % this may cause the models to perform worse because their limited
    % ability for extrapolating across targets
    td1 = td1(randperm(length(td1)));
    trials_p_fold = floor(length(td1)/n_folds);
    
    for f = 1:n_folds
        trials_test = (f-1)*trials_p_fold+1 : f*trials_p_fold;
        trials_train = setdiff(1:length(td1),trials_test);
        td_test = td1(trials_test);
        td_train = td1(trials_train);
        
        [~, mod_info_xval] = getModel(td_train,mod_params);
        % [R2train(f,:), ~] = testModel(td_train,mod_info_xval);
        
        % compute the cross-validated R2
        [R2xval(f,:), ~] = testModel(td_test,mod_info_xval);
    end
    % Save the mean
    withinR2(c,:) = mean(R2xval,1);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now test the model on the data from day 2
        
    [R2diff, td2] = testModel(td2,mod_info);
    acrossR2(c,:) = R2diff;
    
    n_units(c) = size(td1(1).M1_unit_guide,1);
end



% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PREDICTION STABILITY OVER DAYS

fitX = polyfit(diff_days',acrossR2(:,1),1);
fitXy = polyval(fitX,[min(diff_days) max(diff_days)]);
fitY = polyfit(diff_days',acrossR2(:,2),1);
fitYy = polyval(fitY,[min(diff_days) max(diff_days)]);


figure, hold on
plot(diff_days,acrossR2(:,1),'.r','linestyle','none','markersize',30)
plot(diff_days,acrossR2(:,2),'.b','linestyle','none','markersize',30)
ylim([0 1]),legend('X vel','Y vel'),legend boxoff
set(gca,'TickDir','out','FontSize',14)
xlabel('Days since decoder training')
ylabel('Velocity predictions (R^2)')
plot([min(diff_days) max(diff_days)],fitXy,'r','linewidth',2)
plot([min(diff_days) max(diff_days)],fitYy,'b','linewidth',2)
title('Across-day prediction')



% Same plot combining X and Y
fitXY = polyfit(diff_days',mean(acrossR2,2),1);
fitXYy = polyval(fitXY,[min(diff_days) max(diff_days)]);

% add within-day
fitXYw = polyfit(diff_days',mean(withinR2,2),1);
fitXYyw = polyval(fitXYw,[min(diff_days) max(diff_days)]);


figure, hold on
plot(diff_days,mean(acrossR2,2),'.k','linestyle','none','markersize',30)
ylim([0 1]), title('Predictions using threshold crossings')
set(gca,'TickDir','out','FontSize',14)
xlabel('Days since decoder training')
ylabel('Velocity predictions (R^2)')
% adding within day
plot(diff_days,mean(withinR2,2),'.r','linestyle','none','markersize',30)
plot([min(diff_days) max(diff_days)],fitXYy,'k','linewidth',2)
plot([min(diff_days) max(diff_days)],fitXYyw,'r','linewidth',2)
legend('Across-day','Within-day'),legend boxoff



% Plot number of common units
fitUnits = polyfit(diff_days,n_units,1);
fitUnitsy = polyval(fitUnits,[min(diff_days) max(diff_days)]);

figure, hold on
plot(diff_days,n_units,'.','color',[.6 .6 .6],'linestyle','none','markersize',30)
plot([min(diff_days) max(diff_days)],fitUnitsy,'color',[.6 .6 .6],'linewidth',2)
set(gca,'TickDir','out','FontSize',14)
xlabel('Days since decoder training')
ylabel('Number units in decoder')
ylim([0 96])



% Within-day decoders
figure, hold on
plot(diff_days,withinR2(:,1),'.r','linestyle','none','markersize',30)
plot(diff_days,withinR2(:,2),'.b','linestyle','none','markersize',30)
ylim([0 1]),legend('X vel','Y vel'),legend boxoff
set(gca,'TickDir','out','FontSize',14)
xlabel('Days since decoder training')
ylabel('Velocity predictions (R^2)')
% plot([min(diff_days) max(diff_days)],fitXy,'r')
% plot([min(diff_days) max(diff_days)],fitYy,'b')
title('Within-day prediction')