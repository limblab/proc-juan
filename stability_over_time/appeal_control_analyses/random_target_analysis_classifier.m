clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');
folder_name = 'Random Target Analysis v2';

save_figs = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameters

ctr_targ = [-2 2 -32 -27]; % [x low x high y low y high]
ang_width = pi/16; % how wide of a window to look, +/- targ direction

monkey = 'Chewie';
arrays = {'M1','PMd'};

% for classifier do this
idx_start_plan = {'idx_go_cue',-10};
idx_end_plan = {'idx_go_cue',20};
% idx_start_plan = {'idx_movement_on',-20};
% idx_end_plan = {'idx_movement_on',1};

idx_start_move = {'idx_movement_on',-2};
idx_end_move = {'idx_movement_on',11};


avg_dims = 1:4;

pars.align_latent_params.xval_yn        = false;
pars.align_latent_params.n_folds        = 6;
pars.align_latent_params.method         = 'cca'; % 'cca' 'procrustes'
pars.align_latent_params.n_shuff        = 100; % n shuffles for within day ceiling
pars.align_latent_params.save_fig       = false;
pars.align_latent_params.top_lv_plot    = 4;

pars.stab_behav.signal                  = 'vel';
pars.stab_behav.trial_avg               = false;
pars.stab_behav.save_fig                = true;

% Method
pars.class_params.method                = 'Bayes'; % 'NN' 'Bayes'
% Inputs and Outputs
pars.class_params.out                   = 'target_direction';
% History?
pars.class_params.hist_bins             = 0;
% Folds for multi-fold cross-validation
pars.class_params.n_folds               = 100;
pars.class_params.num_test_trials       = 1; % number of trials per target for test set
% A couple other slightly redundant definitions
pars.class_params.idx_start             = 'start';
pars.class_params.idx_end               = 'end';
pars.class_params.idx_start_classify    = 'start';
pars.class_params.idx_end_classify      = 'end';
pars.class_params.save_fig              = true;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep data
% load random target data
load('/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_RT_CS_2016-10-21.mat');
td_rt_store =  trial_data;
% load center out data
load('/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/Chewie/Chewie_CO_CS_2016-10-21.mat');
td_co_store = trial_data;
clear trial_data;

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
td_co_move = trimTD(td,idx_start_move,idx_end_move);
td_co_plan = trimTD(td,idx_start_plan,idx_end_plan);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build fake center out data from random target

td = td_rt;
% compute speed
td = getNorm(td,struct('signals','vel','norm_name','speed'));
% remove incomplete go cues
for trial = 1:length(td)
    bad_idx = isnan(td(trial).idx_go_cue);
    td(trial).idx_go_cue(bad_idx) = [];
    td(trial).target_center(bad_idx,:) = [];
end
% split based on go cue
td = splitTD(td,struct('split_idx_name','idx_go_cue','linked_fields','target_center','extra_bins',[15 15]));
% get movement onset and remove the trials where it doesn't find a good one
td = getMoveOnsetAndPeak(td);
td = removeBadTrials(td,struct('remove_nan_idx',true));

% store this one for later
td_rt_move = trimTD(td,idx_start_move,idx_end_move);


% find all trials where go cue occurs while monkey is in center
good_idx = false(size(td));
for trial = 1:length(td)
    idx = td(trial).idx_movement_on;
    good_idx(trial) = td(trial).pos(idx,1) > ctr_targ(1) & ...
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

rt_trial_idx = good_idx & ~bad_idx;

td = td(rt_trial_idx);
td_rt_plan = trimTD(td,idx_start_plan,idx_end_plan);
td_rt_move = td_rt_move(rt_trial_idx);
% td_co_move = td_co_move(1:length(td_rt_move));


% dim reduce
td_co_plan = dimReduce(td_co_plan,'M1_spikes');
td_co_plan = dimReduce(td_co_plan,'PMd_spikes');
td_rt_plan = dimReduce(td_rt_plan,'M1_spikes');
td_rt_plan = dimReduce(td_rt_plan,'PMd_spikes');

td_co_move = dimReduce(td_co_move,'M1_spikes');
td_co_move = dimReduce(td_co_move,'PMd_spikes');
td_rt_move = dimReduce(td_rt_move,'M1_spikes');
td_rt_move = dimReduce(td_rt_move,'PMd_spikes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build classifier on CO data
array = 'PMd';
switch lower(array)
    case 'm1'
        mani_dims = 1:10;
    case 'pmd'
        mani_dims = 1:16;
end
pars.align_latent_params.mani_dims      = mani_dims;
pars.class_params.mani_dims             = mani_dims;
pars.align_latent_params.signals        = [array '_pca'];
pars.class_params.manifold              = [array '_pca'];

td_co = td_co_plan;
td_rt = td_rt_plan;

master_td = catTDs(td_co,td_rt);
master_td = equalNbrTrialsSessions(master_td);
% res = classify_across_days(master_td,'align',pars.class_params);

dates = unique({master_td.date});
[~,td1] = getTDidx(master_td,'date',dates{1});
[~,td2] = getTDidx(master_td,'date',dates{2});

targs = unique([td2.target_direction]);

f = figure; hold all;
plot(NaN,NaN,'b-','LineWidth',2);
plot(NaN,NaN,'r-','LineWidth',2);

how_many = abs(idx_start_plan{2})+abs(idx_end_plan{2})-5;

td = td2;
p = zeros(length(td),8,how_many);
pred = zeros(length(td),how_many);
for n = 1:how_many
    td_temp = trimTD(td,{'start',n},{'start',n+5});
    X = getSig(binTD(td_temp,'average'),{[array '_pca'],pars.align_latent_params.mani_dims});
    X = zscore(X);
    Y = [td.target_direction]';
    
    clas = fitcnb(X,Y,'ClassNames',targs, ...
        'DistributionNames','normal');
    
    [pred(:,n),p(:,:,n)] = predict(clas,X);
end
perc_correct(1) = 100*sum(pred(:,n)==[td.target_direction]')/size(pred,1);
% now plot dynamics
for trial = 1:length(td)
    if pred(trial,1) ~= td(trial).target_direction && pred(trial,end) == td(trial).target_direction
        idx = targs == td(trial).target_direction;
        plot( (30*(1:how_many)+idx_start_plan{2}*30+180)/1000, ...
            squeeze(p(trial,idx,:)),'r-','LineWidth',1.5);
    end
end

td = td1;
p = zeros(length(td),8,how_many);
pred = zeros(length(td),how_many);
for n = 1:how_many
    td_temp = trimTD(td,{'start',n},{'start',n+5});
    X = getSig(binTD(td_temp,'average'),{[array '_pca'],pars.align_latent_params.mani_dims});
    X = zscore(X);
    Y = [td.target_direction]';
    
    clas = fitcnb(X,Y,'ClassNames',targs, ...
        'DistributionNames','normal');
    
    [pred(:,n),p(:,:,n)] = predict(clas,X);
end
perc_correct(2) = 100*sum(pred(:,n)==[td.target_direction]')/size(pred,1);
% now plot dynamics
for trial = 1:length(td)
    if pred(trial,1) ~= td(trial).target_direction
        idx = targs == td(trial).target_direction;
        plot( (30*(1:how_many)+idx_start_plan{2}*30+180)/1000, ...
            squeeze(p(trial,idx,:)),'b-','LineWidth',1.5);
    end
end

% doctor it up
plot([0 0],[0 1],'k--');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Time (s)');
ylabel('Probability of correct target');
axis tight;
set(gca,'YLim',[0 1]);

h = legend({'CO','RT'},'Location','SouthEast');
set(h,'Box','off');

if save_figs
    fn = [monkey '_COvsRT_ClassifierDynamics'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end

f =  figure;
hold all;
bar(perc_correct([2 1]));
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',[1 2],'XTickLabel',{'CO','RT'});
ylabel('Percent Correct')


if save_figs
    fn = [monkey '_COvsRT_ClassifierPerformance'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end


