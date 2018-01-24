%
%  Params that look good: bin size = 50 ms; SD = 0.1; 2 history bins.
%  Chewie: 8D; Mihili: 
%

clear all, % close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters / definition
monkey          = 'mihili'; % 'mihili'; 'chewie'

unsorted_yn     = false;

spiking_inputs  = {'M1_spikes'}; % 'PMd_spikes'

% The behavioral variable to be predicted
dec_var         = 'vel';
% Decoder inputs
dec_input       = 'aligned_data'; % 'aligned_data'; 'raw_data'

% When we start looking at each trial
idx_start       = {'idx_go_cue',0}; % {'idx_movement_on',0}; % {'idx_go_cue',0}
% When we stop looking at each trial --if empty, the duration of the shortest trial 
idx_end         = {'idx_go_cue',18}; % {''}; % FOR CHEWIE: {'idx_go_cue',10}
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

files_mihili = { ...edit
                'Mihili_CO_FF_2014-02-03.mat', ...
                'Mihili_CO_FF_2014-02-17.mat', ...
                'Mihili_CO_FF_2014-02-18.mat', ...
                'Mihili_CO_VR_2014-03-03.mat', ...
                'Mihili_CO_VR_2014-03-04.mat', ...
                'Mihili_CO_VR_2014-03-06.mat', ...
                'Mihili_CO_FF_2014-03-07.mat' ...
                };
 

% override monkey name and go to path
here        = pwd;
switch monkey
    case 'chewie'
        if ~unsorted_yn
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie/Unsorted')
            for f = 1:length(files_chewie), files_chewie{f} = [files_chewie{f}(1:end-4) '_unsorted.mat']; end
        end
        files = files_chewie;
    case 'mihili'
        if ~unsorted_yn 
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili')
        else
            cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili/Unsorted')
            for f = 1:length(files_mihili), files_mihili{f} = [files_mihili{f}(1:end-4) '_unsorted.mat']; end
        end
        files = files_mihili;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data, unsorting the channels. Baseline trials only

if ~unsorted_yn
    master_td   = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@sqrtTransform,spiking_inputs}, ...
                            {@smoothSignals,struct('signals',spiking_inputs,'calc_fr',true,'kernel_SD',kernel_SD)}, ...
                            {@getPCA,struct('signals',spiking_inputs)}, ...
                            {@binTD,n_bins} ...
                            );
                        
    % some trials have the problem that 
    
else      
    master_td   = loadTDfiles(  files, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@mergeUnits});
                        
%     master_td   = loadTDfiles(  files, ...
%                             {@getTDidx,'epoch','BL'}, ...
%                             {@getTDidx,'result','R'}, ...
%                             {@mergeUnits},...
%                             {@sqrtTransform,spiking_inputs}, ...
%                             {@smoothSignals,struct('signals',spiking_inputs,'calc_fr',true,'kernel_SD',kernel_SD)}, ...
%                             {@getPCA,struct('signals',spiking_inputs)}, ...
%                             {@binTD,n_bins} ...
%                             );

    master_td   = getCommonUnits(master_td);
end

                        
% go back to where you were path-wise
cd(here);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get some meta info about the data

% get the sessions
sessions    = unique({master_td.date});

% get the targets
targets     = unique(cell2mat({master_td.target_direction}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete bad channels if there are any
% --> this is only done when using sorted data; for unsorted, we assume we
% want to do analyses with the threshold crossings, so we leave everything
if ~unsorted_yn
    for s = 1:length(sessions)
        this_s = getTDidx(master_td,{'date',sessions{s}});
        master_td(this_s) = removeBadNeurons(master_td(this_s));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do PCA --only necessary for unsorted datafiles, because in this case we
% don't do PCA during the loading (because we want to have the same
% electrodes in all the sessions)
if unsorted_yn
    
    for s = 1:length(sessions)
        this_s = getTDidx( master_td,{'date',sessions{s}} );
        master_td_new(this_s) = getPCA( master_td(this_s), struct(...
                                    'signals',spiking_inputs,...
                                    'sqrt_transform',true,...
                                    'do_smoothing',true,...
                                    'kernelSD',0.1 ));
    end
    
    master_td = master_td_new;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim the trials

% get the minimum trial duration, if needed for idx_end
if isempty(idx_end{1})
    min_mov_duration = min( arrayfun( @(x) x.idx_trial_end-x.idx_movement_on, master_td) );
    idx_end = {'idx_movement_on',min_mov_duration};
end    

master_td   = trimTD( master_td, idx_start, idx_end );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep the same number of trials for all the targets and sessions --this is
% important to align the neural mode dynamics with CCA

master_td   = equalNbrTrialsSessions(master_td);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do CCA of the neural trajectories across sessions, to align the dynamics
% and see how similar they are

% get all pairwise comparisons of sessions
comb_sessions   = nchoosek(1:length(sessions),2);

% do
disp('Aligning neural mode dynamics across days')
for c = 1:size(comb_sessions,1)
    
    % get the data
    trials1 = getTDidx( master_td, 'date', sessions{comb_sessions(c,1)} );
    trials2 = getTDidx( master_td, 'date', sessions{comb_sessions(c,2)} );
    
    %trials2 = trials2(randperm(length(trials2)));

    switch spiking_inputs{1}
        case 'M1_spikes'
            manifold = 'M1_pca';
        case 'PMd_spikes'
            manifold = 'PMd_pca';
    end
    % compare dynamics with CCA
    cca_info(c) = compDynamics( master_td, manifold, trials1, trials2, mani_dims );
    
    % compare dynamics with good old forrelations
    corr_info(c) = corrDynamics( master_td, manifold, trials1, trials2, mani_dims );
end


% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TEST STABILITY OVER TIME OF DECODERS THAT USE LATENT SIGNALS AS INPUTS:
% BUILD A DECODER ON ONE DAY AND TEST IT ON ANOTHER. Do for all pairs of
% sessions


disp('Testing decoder stability');
disp(['Decoding: ' dec_var]);

for c = 1:size(comb_sessions,1)

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add aligned neural mode dynamics to trial_data
    
    % get two TD structs, one per session that we are comparing. The code
    % will build a model on TD1 and test it on TD2
    [trials1, td1] = getTDidx(master_td,'date',sessions{comb_sessions(c,1)});
    [trials2, td2] = getTDidx(master_td,'date',sessions{comb_sessions(c,2)});
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SANTIY CHECK: the CC between the aligned & unaligned data is all 1 as
    % it should
    %
    x = get_vars(td1,{manifold,mani_dims});
    y = get_vars(td2,{manifold,mani_dims});
    u = cca_info(c).U;
    v = cca_info(c).V;
    [~,~,cc_xu] = canoncorr(x,u);
    [~,~,cc_yv] = canoncorr(y,v);
    if min(cc_xu) < .99, error(['CC between raw and aligned data ~1 for c = ' num2str(c)]); end
    if min(cc_yv) < .99, error(['CC between raw and aligned data ~1 for c = ' num2str(c)]); end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Separate concatenated neural mode dynamics into trials
    bins_p_trial = size(td1(1).pos,1);
    
    for t = 1:length(trials1)
        be = bins_p_trial*(t-1)+1;
        en = bins_p_trial*t;        
        td1(t).(['align_' manifold]) = cca_info(c).U(be:en,:);
        td2(t).(['align_' manifold]) = cca_info(c).V(be:en,:);
    end
    
    
    % duplicate and shift --to have history in the model
    td1 = dupeAndShift(td1,{['align_' manifold],hist_bins});
    td2 = dupeAndShift(td2,{['align_' manifold],hist_bins});

    
    % build linear model
    mod_params.model_type = 'linmodel';
    mod_params.out_signals = dec_var;
    switch dec_input
        case 'aligned_data'
            mod_params.in_signals = {['align_' manifold '_shift'],1:length(mani_dims)*hist_bins};
        case 'raw_data'
            td1 = dupeAndShift(td1,{manifold,hist_bins});
            td2 = dupeAndShift(td2,{manifold,hist_bins});
            mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
        otherwise
            error('wrong decoder input');
    end
    [td1, mod_info] = getModel(td1,mod_params);

    
    % Do multi-fold cross-validation of the within-day predictions
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
        
    [R2diff td2] = testModel(td2,mod_info);
    acrossR2(c,:) = R2diff;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % There seems to be a gain offset -> compute it!
    
    Yacross = get_vars(td2,{dec_var,1:size(td2(1).(dec_var),2)});
    Yhat_across = get_vars(td2,{'linmodel_default',1:size(td2(1).(dec_var),2)});
    
    gains(c,:) = computeGain(Yacross,Yhat_across);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % As control, compute how well a decoder generalizes if we don't align
    % the data
    
    % duplicate and shift --to have history in the model
    ctrl_td1 = dupeAndShift(td1,{manifold,hist_bins});
    ctrl_td2 = dupeAndShift(td2,{manifold,hist_bins});

    % build linear model
    ctrl_mod_params.model_type = 'linmodel';
    ctrl_mod_params.out_signals = dec_var;
    ctrl_mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
    
    [ctrl_td1, ctrl_mod_info] = getModel(ctrl_td1,ctrl_mod_params);
    
    % use the model from day 1 to predict day 2
    [R2ctrl ctrl_td2] = testModel(ctrl_td2,ctrl_mod_info);
    ctrlR2(c,:) = R2ctrl;
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAME ANALYSIS USING SPIKES (SINGLE NEURONS)


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
    withinR2_spikes(c,:) = mean(R2xval,1);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now test the model on the data from day 2
        
    [R2diff, td2] = testModel(td2,mod_info);
    acrossR2_spikes(c,:) = R2diff;
    
    n_units(c) = size(td1(1).M1_unit_guide,1);
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SANITY CHECKS

% -------------------------------------------------------------------------
% BUILD A MODEL ON THE FIRST HALF OF EACH SESSION, AND TEST IT ON THE OTHER
% HALF --Working

% for s = 1:length(sessions)
%     
%     % retrieve trials for this session
%     [~, this_td] = getTDidx(master_td,{'date',sessions{s}});
%     
%     % duplicate and shift to have history
%     this_td = dupeAndShift(this_td,{manifold,hist_bins});
%     
%     % use same model params 
%     this_mod_params.model_type = mod_params.model_type;
%     this_mod_params.out_signals = mod_params.out_signals;
%     this_mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
%     
%     % split trials in two randomly -necessary because the trials to the
%     % same target are ordered
%     random_trials = randperm(length(this_td));
%     
%     this_train_trials = random_trials(1:length(random_trials)/2);
%     this_test_trials = random_trials(length(random_trials)/2+1:end);
%     
%     this_train_td = this_td(this_train_trials);
%     this_test_td = this_td(this_test_trials);
%     
%     % train model on first half on the data
%     [this_train_td, this_mod_info] = getModel(this_train_td,this_mod_params);
% 
%     % and get the within-fit R2
%     this_Y = get_vars(this_train_td,{dec_var,1:size(this_train_td(1).(dec_var),2)});
%     this_Yhat = get_vars(this_train_td,{'linmodel_default',1:size(this_train_td(1).(dec_var),2)});
%    
%     this_within_R2 = compute_r2(this_Y,this_Yhat);
%     
%     sanity_within_R2(s,:) = this_within_R2;
%     
%     % test model on second half of data
%     [this_across_R2 this_test_td] = testModel(this_test_td,this_mod_info);
%     sanity_across_R2(s,:) = this_across_R2;
% end
% 
% 
% figure,hold on
% plot(sanity_across_R2(:,1),'b','linewidth',2); 
% plot(sanity_across_R2(:,2),'k','linewidth',2); 
% plot(sanity_within_R2(:,1),'c','linewidth',2); 
% plot(sanity_within_R2(:,2),'color',[.6 .6 .6],'linewidth',2); 
% set(gca,'TickDir','out','FontSize',14), box off
% ylim([0 1]),xlim([0 size(sanity_across_R2,1)+1])
% legend('X train','Y train','X test','Y test','Location','Southwest'),legend boxoff
% xlabel('Session'); ylabel('Model accuracy (R^2)')


% -------------------------------------------------------------------------
% BUILD A MODEL USING THRESHOLD CROSSINGS AND USE IT ACROSS DAYS






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within-day predictions using latent signals

figure,hold on, 
plot(unique(round(withinR2(:,1)*10000)/10000,'rows'),'.b','markersize',32,'linewidth',2)
plot(unique(round(withinR2(:,2)*10000)/10000,'rows'),'.r','markersize',32,'linewidth',2)
set(gca,'TickDir','out','FontSize',14), box off
legend('X','Y','Location','SouthEast'),legend('boxoff')
ylabel('Model R^2'), ylim([0 1])
set(gca,'XTick',1:length(sessions)-1,'XTickLabel',sessions(1:end-1),'XTickLabelRotation',45)
xlim([0 length(sessions)])
switch dec_input
    case 'aligned_data'
        title(['Within-day predictions from an aligned manifold with D = ' num2str(length(mani_dims))]);
    case 'raw_data'
        title(['Within-day predictions from a manifold with D = ' num2str(length(mani_dims))]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Across day predictions - separated per predicted variable

% % Linear fits for plots
% for d = 1:size(acrossR2,2)
%     lin_fit(d,:) = polyfit(diff_days',acrossR2(:,d),1);
%     coordplot(d,:) = polyval(lin_fit(d,:),[min(diff_days)-1, max(diff_days)+1]);
%     lin_fit_ctrl(d,:) = polyfit(diff_days',ctrlR2(:,d),1);
%     coordplot_ctrl(d,:) = polyval(lin_fit_ctrl(d,:),[min(diff_days)-1, max(diff_days)+1]);
% end
%    
% figure,hold on, 
% % unaligned data
% plot(diff_days,ctrlR2(:,1),'o','color',[.6 .6 .6],'markersize',5,'linewidth',2)
% plot(diff_days,ctrlR2(:,2),'.','color',[.6 .6 .6],'markersize',20,'linewidth',2)
% % aligned data
% plot(diff_days,acrossR2(:,1),'.b','markersize',32,'linewidth',2)
% plot(diff_days,acrossR2(:,2),'.r','markersize',32,'linewidth',2)
% plot([min(diff_days),max(diff_days)],coordplot(1,:),'b','linewidth',1)
% plot([min(diff_days),max(diff_days)],coordplot(2,:),'r','linewidth',1)
% plot([min(diff_days),max(diff_days)],coordplot_ctrl(1,:),'color',[.6 .6 .6],'linestyle','-.','linewidth',1)
% plot([min(diff_days),max(diff_days)],coordplot_ctrl(2,:),'color',[.6 .6 .6],'linewidth',1)
% set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+10])
% legend('X_{unalinged}','Y_{unaligned}','X','Y','Location','NorthEast'),legend('boxoff'),title(monkey)
% ylabel('Model stability over sessions (R^2)'),xlabel('Days between sessions'),ylim([0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Across day predictions - one data point per day (mean across predicted
% variables)

% linear fits
lin_fit_pool = polyfit(diff_days',mean(acrossR2,2),1);
coordplot_pool = polyval(lin_fit_pool,[min(diff_days)-1, max(diff_days)+1]);
lin_fit_ctrl_pool = polyfit(diff_days',mean(ctrlR2,2),1);
coordplot_ctrl_pool = polyval(lin_fit_ctrl_pool,[min(diff_days)-1, max(diff_days)+1]);

figure,hold on, 
% unaligned data
plot(diff_days,mean(ctrlR2,2),'.','color',[.6 .6 .6],'markersize',24,'linewidth',2)
% aligned data
plot(diff_days,mean(acrossR2,2),'.b','markersize',32,'linewidth',2)
plot([min(diff_days),max(diff_days)],coordplot_pool,'b','linewidth',1)
plot([min(diff_days),max(diff_days)],coordplot_ctrl_pool,'color',[.6 .6 .6],'linestyle','-.','linewidth',1)
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+1])
legend('unaligned','aligned','Location','NorthEast'),legend('boxoff'),title(monkey)
ylabel('Velocity predicitons based on latent signals (R^2)'),xlabel('Days between sessions'),ylim([0 1])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Across day predictions + within-day predictions - one data point per day
% (mean across predicted variables)

% linear fits --copied and pasted from above
lin_fit_pool = polyfit(diff_days',mean(acrossR2,2),1);
coordplot_pool = polyval(lin_fit_pool,[min(diff_days)-1, max(diff_days)+1]);
lin_fit_ctrl_pool = polyfit(diff_days',mean(ctrlR2,2),1);
coordplot_ctrl_pool = polyval(lin_fit_ctrl_pool,[min(diff_days)-1, max(diff_days)+1]);
lin_fit_within_pool = polyfit(diff_days',mean(withinR2,2),1);
coord_within_pool = polyval(lin_fit_within_pool,[min(diff_days)-1, max(diff_days)+1]);

figure,hold on, 
% unaligned data
plot(diff_days,mean(ctrlR2,2),'.','color',[.6 .6 .6],'markersize',32,'linewidth',2)
% aligned data
plot(diff_days,mean(acrossR2,2),'.b','markersize',32,'linewidth',2)
% within-day data
plot(diff_days,mean(withinR2,2),'.k','markersize',32,'linewidth',2)
plot([min(diff_days),max(diff_days)],coordplot_pool,'b','linewidth',1)
plot([min(diff_days),max(diff_days)],coordplot_ctrl_pool,'color',[.6 .6 .6],'linewidth',1)
plot([min(diff_days),max(diff_days)],coord_within_pool,'k','linewidth',1)
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+1])
legend('unaligned','aligned','wthin-day','Location','NorthEast'),legend('boxoff'),title(monkey)
ylabel('Velocity predicitons based on latent signals (R^2)'),xlabel('Days between sessions'),ylim([0 1])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within-day vs Across-day predictions with decoders that use spikes as
% inputs - one data point per day - everything xval


% Same plot combining X and Y
fitXY = polyfit(diff_days',mean(acrossR2_spikes,2),1);
fitXYy = polyval(fitXY,[min(diff_days) max(diff_days)]);

% add within-day
fitXYw = polyfit(diff_days',mean(withinR2_spikes,2),1);
fitXYyw = polyval(fitXYw,[min(diff_days) max(diff_days)]);


figure, hold on
plot(diff_days,mean(acrossR2_spikes,2),'.k','linestyle','none','markersize',30)
ylim([0 1]), title('Predictions using threshold crossings')
set(gca,'TickDir','out','FontSize',14)
xlabel('Days since decoder training')
ylabel('Velocity predictions based on single neurons (R^2)')
% adding within day
plot(diff_days,mean(withinR2_spikes,2),'.r','linestyle','none','markersize',30)
plot([min(diff_days) max(diff_days)],fitXYy,'k','linewidth',2)
plot([min(diff_days) max(diff_days)],fitXYyw,'r','linewidth',2)
legend('Across-day','Within-day'),legend boxoff




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within-day vs Across-day predictions with decoders that use latent
% signals and spikes as inputs - one data point per day - everything xval


figure,hold on, 
% % unaligned, latent signals
% plot(diff_days,mean(ctrlR2,2),'.','color',[.6 .6 .6],'markersize',32,'linewidth',2)
% across days: aligned, latent signals
plot(diff_days,mean(acrossR2,2),'.b','markersize',32,'linewidth',2)
% within-day, latent signals
plot(diff_days,mean(withinR2,2),'.k','markersize',32,'linewidth',2)
% across days, spikes
plot(diff_days,mean(acrossR2_spikes,2),'.r','markersize',32,'linewidth',2)
% within-day, spikes
plot(diff_days,mean(withinR2_spikes,2),'color',[.6 .6 .6],'marker','.','linestyle','none','markersize',32)
plot([min(diff_days),max(diff_days)],coordplot_pool,'b','linewidth',1)
% plot([min(diff_days),max(diff_days)],coordplot_ctrl_pool,'color',[.6 .6 .6],'linewidth',1)
plot([min(diff_days),max(diff_days)],coord_within_pool,'k','linewidth',1)
plot([min(diff_days) max(diff_days)],fitXYy,'r','linewidth',2)
plot([min(diff_days) max(diff_days)],fitXYyw,'color',[.6 .6 .6],'linewidth',2)
legend('Across-day latent signals','Within-day latent signals','Across-day neurons','Within-day neurons',...
    'Location','SouthWest'),
legend boxoff
title(monkey)
set(gca,'TickDir','out','FontSize',14)
ylim([0 1]), ylabel('Velocity predictions (R^2)')
xlabel('Days since decoder training')
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+1])
legend('across latent','within-day latent','across neurons','within-day neurons','Location','SouthWest'),
legend('boxoff'),title(monkey)
ylabel('Velocity predicitons (R^2)'),
xlabel('Days between sessions'),ylim([0 1])





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Change in gain as a function of time, to show that something in the
% % mapping changes
% 
% figure,hold on
% plot(diff_days,gains(:,1),'.b','markersize',32)
% plot(diff_days,gains(:,2),'.r','markersize',32)
% set(gca,'TickDir','out','FontSize',14), box off
% xlabel('Days between sessions')
% ylabel('Magnitude of actual / predicted signals'),ylim([0 1])
% legend('X','Y','Location','NorthEast'),legend('boxoff'),title(monkey)
% title(monkey)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Chain in gain vs model accuracy
% 
% figure,hold on
% plot(acrossR2(:,1),gains(:,1),'.b','markersize',32)
% plot(acrossR2(:,2),gains(:,2),'.r','markersize',32)
% set(gca,'TickDir','out','FontSize',14), box off
% xlabel('Model stability over sessions (R^2)'),xlim([0 1])
% ylabel('Magnitude of actual / predicted signals'),ylim([0 1])
% legend('X','Y','Location','NorthWest'),legend('boxoff'),title(monkey)
