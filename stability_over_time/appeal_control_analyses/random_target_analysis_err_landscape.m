clear all;
close all;
clc;

% params
array = 'M1';

switch lower(array)
    case 'pmd'
        idx_start = {'idx_movement_on',-13};
        idx_end = {'idx_movement_on',0};
        pars.align_latent_params.mani_dims      = 1:16;
    case 'm1'
        idx_start = {'idx_movement_on',-2};
        idx_end = {'idx_movement_on',11};
        pars.align_latent_params.mani_dims      = 1:10;
end

idx_start_move = {'idx_movement_on',0};
idx_end_move = {'idx_movement_on',13};

avg_dims = 1:4;

pars.align_latent_params.xval_yn        = false;
pars.align_latent_params.n_folds        = 6;
pars.align_latent_params.method         = 'cca'; % 'cca' 'procrustes'
pars.align_latent_params.n_shuff        = 100; % n shuffles for within day ceiling
pars.align_latent_params.signals        = [array '_pca'];

pars.align_latent_params.save_fig       = false;
pars.align_latent_params.top_lv_plot    = ceil(length(pars.align_latent_params.mani_dims)/3);

pars.stab_behav.signal                  = 'vel';
pars.stab_behav.trial_avg               = false;
pars.stab_behav.save_fig                = true;



%% prep data
% load random target data
load('/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_RT_CS_2016-10-21.mat');
td_rt_store =  trial_data;
% load center out data
load('/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_CO_CS_2016-10-21.mat');
td_co_store = trial_data;
clear trial_data;

%% go for it
% define parameters
clear run_vals;
run_vals(1).ctr_targ = [-2 2 -32 -28];  % [x low x high y low y high]
run_vals(1).ang_width = pi/16;
run_vals(2).ctr_targ = [-2 2 -32 -28];  % [x low x high y low y high]
run_vals(2).ang_width = pi/6;
run_vals(3).ctr_targ = [-2 2 -32 -28];  % [x low x high y low y high]
run_vals(3).ang_width = pi/4;
run_vals(4).ctr_targ = [-4 4 -34 -26];  % [x low x high y low y high]
run_vals(4).ang_width = pi/4;
run_vals(5).ctr_targ = [-8 8 -38 -22];  % [x low x high y low y high]
run_vals(5).ang_width = pi/4;
run_vals(6).ctr_targ = [-10 10 -40 -20];  % [x low x high y low y high]
run_vals(6).ang_width = pi/4;

results = zeros(length(run_vals),4);
err = zeros(1,length(run_vals));
for iRun = 1:length(run_vals)
    disp(iRun);
    
    ctr_targ = run_vals(iRun).ctr_targ;
    ang_width = run_vals(iRun).ang_width;
    
    td = td_rt_store;
    td = binTD(td,3);
    td = sqrtTransform(td);
    td =  smoothSignals(td,struct('signals','M1_spikes','width',0.05,'calc_fr',true));
    td =  smoothSignals(td,struct('signals','PMd_spikes','width',0.05,'calc_fr',true));
    td_rt =  td;
    
    td = td_co_store;
    td = binTD(td,3);
    td = sqrtTransform(td);
    td =  smoothSignals(td,struct('signals','M1_spikes','width',0.05,'calc_fr',true));
    td =  smoothSignals(td,struct('signals','PMd_spikes','width',0.05,'calc_fr',true));
            td = removeBadTrials(td,struct('remove_nan_idx',false,  ...
        'ranges',{{idx_end_move{1},'end',[idx_end_move{2}+1 Inf]}}));
    td_co_move = td;
    td = trimTD(td,idx_start,idx_end);
    td_co =  td;
    
    td_co_move = trimTD(td_co_move,idx_start_move,idx_end_move);
    
    % build fake center out data from random target
    
    td = td_rt;
    td = getNorm(td,struct('signals','vel','norm_name','speed'));
    % remove incomplete go cues
    for trial = 1:length(td)
        bad_idx = isnan(td(trial).idx_go_cue);
        td(trial).idx_go_cue(bad_idx) = [];
        td(trial).target_center(bad_idx,:) = [];
    end
    % split based on go cue
    td = splitTD(td,struct('split_idx_name','idx_go_cue','linked_fields','target_center','extra_bins',[20 20]));
    % get movement onset and remove the trials where it doesn't find a good one
    td = getMoveOnsetAndPeak(td);
    td = removeBadTrials(td,struct('remove_nan_idx',true, ...
        'ranges',{{idx_end_move{1},'end',[idx_end_move{2}+1 Inf]}}));
    
    
    % store this one for later
    td_rt_move = td;
    td_rt_move = trimTD(td_rt_move,idx_start_move,idx_end_move);
    
    
    % find all trials where go cue occurs while monkey is in center
    trial_idx = false(size(td));
    for trial = 1:length(td)
        idx = td(trial).idx_movement_on;
        trial_idx(trial) = td(trial).pos(idx,1) > ctr_targ(1) & ...
            td(trial).pos(idx,1) < ctr_targ(2) & ...
            td(trial).pos(idx,2) > ctr_targ(3) & ...
            td(trial).pos(idx,2) < ctr_targ(4);
    end
    
    % now subselect to  only pick ones with the proper angle
    targs = unique([td_co_store.target_direction]);
    bad_idx = false(size(td));
    for trial = 1:length(td)
        good_trial = false(1,length(targs));
        
        % get the movement angle
        th = atan2(td_rt_move(trial).pos(end,2) - td_rt_move(trial).pos(1,2), ...
            td_rt_move(trial).pos(end,1) - td_rt_move(trial).pos(1,1));
        
        for u = 1:length(targs)
            % get angle window
            win = targs(u) + ang_width*[-1 1];
            good_trial(u) = ~isempty(range_intersection(angleDiff(win,th),[0 0]));
        end
        bad_idx(trial) = ~any(good_trial);
        
        temp = bin_angles(th,pi/4);
        if temp == -pi, temp =  pi; end
        td(trial).target_direction = temp;
        td_rt_move(trial).target_direction = temp;
        td(trial).date = datestr(datenum(td(trial).date,'mm-dd-yyyy')+1,'mm-dd-yyyy');
        td_rt_move(trial).date = datestr(datenum(td_rt_move(trial).date,'mm-dd-yyyy')+1,'mm-dd-yyyy');
    end
    
    td = td(trial_idx & ~bad_idx);
    % td([td.idx_movement_on] <= abs(idx_start{2})) = [];
    
    td = trimTD(td,idx_start,idx_end);
    
    td_rt = td;
    
    td_rt_move = td_rt_move(trial_idx & ~bad_idx);
    
%     td_co_move = td_co_move(1:length(td_rt_move));
    
    % dim reduce
    td_co = dimReduce(td_co,'M1_spikes');
    td_co = dimReduce(td_co,'PMd_spikes');
    td_rt = dimReduce(td_rt,'M1_spikes');
    td_rt = dimReduce(td_rt,'PMd_spikes');
    
    % now start doing analyses
    td = catTDs(td_co,td_rt);
    td = equalNbrTrialsSessions(td);
    within_align = align_latent_activity_within_day(td, pars.align_latent_params );
    results(iRun,1) = mean(within_align(1).cc_m(avg_dims));
    results(iRun,2) = mean(within_align(2).cc_m(avg_dims));
    align_results = align_latent_activity(td,pars.align_latent_params);
    results(iRun,3) = mean(align_results.cc(avg_dims));
    
    
    % make the master_td
    master_td = catTDs(td_co_move,td_rt_move);
    master_td = equalNbrTrialsSessions(master_td);
    corr_kin = comp_behavior( master_td, pars.stab_behav );
    results(iRun,4) = mean(corr_kin.r);
    
    dates = unique({master_td.date});
    [~,td1] =  getTDidx(master_td,'date',dates{1});
    [~,td2] =  getTDidx(master_td,'date',dates{2});
    
    err(iRun) = mean(sqrt(mean((getSig(td1,'pos')-getSig(td2,'pos')).^2)));
end

%%
figure;
plot(err, results(:,2:3),'LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'YLim',[0 1]);
ylabel('Canon Corr');
xlabel('"Error" between CO and RT');
legend({'Within RT','Across RT to CO'});

title(array)

%% %%%%%%%
figure('Position',[100 100 500 800]);

td = td_rt_store;

subplot(3,2,1);
hold all;
plot(getSig(td,{'pos',1}),getSig(td,{'pos',2}))
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-11  11], 'YLim',[-41, -19]);
xlabel('X Pos');
ylabel('Y Pos');
title('Random Target (All Movements)');

subplot(3,2,2);
hold all;
plot(getSig(td,{'vel',1}),getSig(td,{'vel',2}))
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-50 50],'YLim',[-50 50]);
xlabel('X Vel');
ylabel('Y Vel');



td = td_rt_move;

subplot(3,2,3);
hold all;
plot(getSig(td,{'pos',1}),getSig(td,{'pos',2}),'.','MarkerSize',1)
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-11  11], 'YLim',[-41, -19]);
xlabel('X Pos');
ylabel('Y Pos');
title('Random Target (Center movements)');

subplot(3,2,4);
hold all;
plot(getSig(td,{'vel',1}),getSig(td,{'vel',2}),'.','MarkerSize',1)
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-50 50],'YLim',[-50 50]);
xlabel('X Vel');
ylabel('Y Vel');


td =  td_co_move;

subplot(3,2,5);
hold all;
plot(getSig(td,{'pos',1}),getSig(td,{'pos',2}),'.','MarkerSize',1)
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-11  11], 'YLim',[-41, -19]);
xlabel('X Pos');
ylabel('Y Pos');
title('Center out');

subplot(3,2,6);
hold all;
plot(getSig(td,{'vel',1}),getSig(td,{'vel',2}),'.','MarkerSize',1)
set(gca,'Box','off','TickDir','out','FontSize',14);
axis square;
set(gca,'XLim',[-50 50],'YLim',[-50 50]);
xlabel('X Vel');
ylabel('Y Vel');




%%
figure;
hold all;
m = mean(results);
s =  std(results);
[hbar,herr] = barwitherr(s,m);
set(hbar,'FaceColor','k','FaceAlpha',0.1,'LineWidth',2)
set(herr,'LineWidth',2)

% do statistics
s = prctile(results(:,3)-results(:,1),[0.5 99.5]);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([1 3],0.98*[1 1],'k','LineWidth',2);
    text(2,0.98,'*','FontSize',30);
end
s = prctile(results(:,2)-results(:,1),[0.5 99.5]);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([1 2],0.94*[1 1],'k','LineWidth',2);
    text(1.5,0.94,'*','FontSize',30);
end
s = prctile(results(:,3)-results(:,2),[0.5 99.5]);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([2 3],0.92*[1 1],'k','LineWidth',2);
    text(2.5,0.92,'*','FontSize',30);
end

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',1:3,'XTickLabel',{'Within CO','Within RT','CO to RT'});
set(gca,'YLim',[0 1]);
s
ylabel('Canon. Corr.');
title(array);


