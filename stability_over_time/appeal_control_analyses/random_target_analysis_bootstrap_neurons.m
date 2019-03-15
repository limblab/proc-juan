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
num_iters = 1000;

ctr_targ = [-2 2 -32 -27]; % [x low x high y low y high]
ang_width = pi/8; % how wide of a window to look, +/- targ direction
min_distance = 0;

monkey = 'Chewie';
arrays = {'M1','PMd'};

idx_start_plan = {'idx_movement_on',-13};
idx_end_plan = {'idx_movement_on',0};

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
    
    % get reach distance
    d = sqrt( ...
        (td_rt_move(trial).pos(end,2) - td_rt_move(trial).pos(1,2)).^2 + ...
        (td_rt_move(trial).pos(end,1) - td_rt_move(trial).pos(1,1)).^2);
    
    bad_idx(trial) = bad_idx(trial) | d < min_distance;
    
    % get the movement angle
    th = atan2(td_rt_move(trial).pos(end,2) - td_rt_move(trial).pos(1,2), ...
        td_rt_move(trial).pos(end,1) - td_rt_move(trial).pos(1,1));
    
    for u = 1:length(targs)
        % get angle window
        win = targs(u) + ang_width*[-1 1];
        good_trial(u) = ~isempty(range_intersection(angleDiff(win,th),[0 0]));
    end
    bad_idx(trial) = bad_idx(trial) | ~any(good_trial);
    
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

% 
% % dim reduce
% td_co_plan = dimReduce(td_co_plan,'M1_spikes');
% td_co_plan = dimReduce(td_co_plan,'PMd_spikes');
% td_rt_plan = dimReduce(td_rt_plan,'M1_spikes');
% td_rt_plan = dimReduce(td_rt_plan,'PMd_spikes');
% 
% td_co_move = dimReduce(td_co_move,'M1_spikes');
% td_co_move = dimReduce(td_co_move,'PMd_spikes');
% td_rt_move = dimReduce(td_rt_move,'M1_spikes');
% td_rt_move = dimReduce(td_rt_move,'PMd_spikes');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now do alignment analyses

array_results = struct();
for iArray = 1:length(arrays)
    array = arrays{iArray};
    
    switch lower(array)
        case 'm1'
            mani_dims = 1:10;
            td_co = td_co_move;
            td_rt = td_rt_move;
        case 'pmd'
            mani_dims = 1:16;
            td_co = td_co_plan;
            td_rt = td_rt_plan;
    end
    pars.align_latent_params.mani_dims      = mani_dims;
    pars.class_params.mani_dims             = mani_dims;
    pars.align_latent_params.signals        = [array '_pca'];
    pars.class_params.manifold              = [array '_pca'];
    
    results = zeros(num_iters,3);
    for iBoot = 1:num_iters
        disp(iBoot)
        td = [];
        % make sure our random bootstrap sample gets enough from each targ
        while length(td) < 5*2*8
            td1 =  td_co(randi(length(td_co),[1,length(td_co)]));
            td2 =  td_rt(randi(length(td_rt),[1,length(td_rt)]));
            
            % random neurons
            idx = randperm(size(td1(1).([array '_spikes']),2));
            for trial = 1:length(td1)
                temp = td1(trial).([array '_spikes']);
                td1(trial).([array '_spikes']) = temp(:,idx(1:floor(length(idx)/2)));
            end
            for trial = 1:length(td2)
                temp = td2(trial).([array '_spikes']);
                td2(trial).([array '_spikes']) = temp(:,idx(floor(length(idx)/2)+1:end));
            end
            
            % dim reduce
            td1 = dimReduce(td1,[array '_spikes']);
            td2 = dimReduce(td2,[array '_spikes']);
            
            td = catTDs(td1,td2);
            td = equalNbrTrialsSessions(td);
        end
        
        within_align = align_latent_activity_within_day(td, pars.align_latent_params );
        results(iBoot,1) = mean(within_align(1).cc_m(avg_dims));
        results(iBoot,2) = mean(within_align(2).cc_m(avg_dims));
        align_results = align_latent_activity(td,pars.align_latent_params);
        results(iBoot,3) = mean(align_results.cc(avg_dims));
    end
    
    array_results.(array) = results;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do behavior correlation for movement
master_td = catTDs(td_co_move,td_rt_move);
master_td = equalNbrTrialsSessions(master_td);
corr_kin        = comp_behavior( master_td, pars.stab_behav );
align_results = align_latent_activity(master_td,pars.align_latent_params);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% latents
if 0
    f = figure;
    
    which_targ = 1:8;
    which_dim = 2;
    
    td_co = td_co_plan;
    td_rt = td_rt_plan;
    for iArray = 1:length(arrays)
        array = arrays{iArray};
        
        
        switch lower(array)
            case 'm1'
                mani_dims = 1:10;
            case 'pmd'
                mani_dims = 1:16;
        end
        pars.align_latent_params.mani_dims      = mani_dims;
        pars.align_latent_params.signals        = [array '_pca'];
        
        master_td = catTDs(td_co,td_rt);
        master_td = equalNbrTrialsSessions(master_td);
        align = align_latent_activity(master_td,pars.align_latent_params);
        
        count = 0;
        for trial = 1:floor(length(master_td)/2)
            master_td(trial).([array '_align']) = align.aligned_info.U(count+1:count + size(master_td(trial).([array '_pca']),1),:);
            count = count + size(master_td(trial).([array '_pca']),1);
        end
        count = 0;
        for trial = floor(length(master_td)/2)+1:length(master_td)
            master_td(trial).([array '_align']) = align.aligned_info.V(count+1:count + size(master_td(trial).([array '_pca']),1),:);
            count = count + size(master_td(trial).([array '_pca']),1);
        end
        dates = unique({master_td.date});
        [~,td_co] = getTDidx(master_td,'date',dates{1});
        [~,td_rt] = getTDidx(master_td,'date',dates{2});
    end
    
    subplot(1,2,1); hold all;
    idx  = getTDidx(td_co,'target_direction',targs(which_targ));
    plot(squeeze(getSigByTrial(td_co(idx),{'PMd_align',which_dim })),'b')
    idx  = getTDidx(td_rt,'target_direction',targs(which_targ));
    plot(squeeze(getSigByTrial(td_rt(idx),{'PMd_align',which_dim })),'r')
    axis tight;
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('Aligned Mode 1');
    title('PMd during Planning')
    
    
    
    td_co = td_co_move;
    td_rt = td_rt_move;
    for iArray = 1:length(arrays)
        array = arrays{iArray};
        
        switch lower(array)
            case 'm1'
                mani_dims = 1:10;
            case 'pmd'
                mani_dims = 1:16;
        end
        pars.align_latent_params.mani_dims      = mani_dims;
        pars.align_latent_params.signals        = [array '_pca'];
        
        master_td = catTDs(td_co,td_rt);
        master_td = equalNbrTrialsSessions(master_td);
        align = align_latent_activity(master_td,pars.align_latent_params);
        
        count = 0;
        for trial = 1:floor(length(master_td)/2)
            master_td(trial).([array '_align']) = align.aligned_info.U(count+1:count + size(master_td(trial).([array '_pca']),1),:);
            count = count + size(master_td(trial).([array '_pca']),1);
        end
        count = 0;
        for trial = floor(length(master_td)/2)+1:length(master_td)
            master_td(trial).([array '_align']) = align.aligned_info.V(count+1:count + size(master_td(trial).([array '_pca']),1),:);
            count = count + size(master_td(trial).([array '_pca']),1);
        end
        dates = unique({master_td.date});
        [~,td_co] = getTDidx(master_td,'date',dates{1});
        [~,td_rt] = getTDidx(master_td,'date',dates{2});
    end
    
    subplot(1,2,2); hold all;
    idx  = getTDidx(td_co,'target_direction',targs(which_targ));
    plot(squeeze(getSigByTrial(td_co(idx),{'M1_align',which_dim })),'b')
    idx  = getTDidx(td_rt,'target_direction',targs(which_targ));
    plot(squeeze(getSigByTrial(td_rt(idx),{'M1_align',which_dim })),'r')
    axis tight;
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('Aligned Mode 1');
    title('M1 during movement')
    
    if save_figs
        fn = [monkey '_COvsRT_LatentDynamics'];
        saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot position traces
f = figure('Position',[100 100 500 800]);

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


if save_figs
    fn = [monkey '_COvsRT_PositionTraces'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot "clouds"
tds = {td_co_move, td_rt_move};
td_names = {'CO','RT'};
epochs = { ...
    %     {'start',0}, {'start',5}; ...
    {'idx_movement_on',-3}, {'idx_movement_on',1}; ...
    {'idx_movement_on',0}, {'idx_movement_on',10}; ...
    };
epoch_names = {'Planning','Movement'};
for iArray = 1:length(arrays)
    f = figure('Position',[100 100 600 600]);
    
    array = arrays{iArray};
    
    ax = [];
    for i = 1:length(tds)
        for j = 1:size(epochs,1)
            td = tds{i};
            td = trimTD(td,epochs{j,1},epochs{j,2});
            td = binTD(td,'average');
            
            targs = unique([td.target_direction]);
            plot_colors = distinguishable_colors(length(targs));
            
            ax = [ax, subplot(size(epochs,1),length(tds),length(tds)*(j-1)+i)]; hold all;
            for trial = 1:length(td)
                targ = targs == td(trial).target_direction;
                
                temp = td(trial).([array '_pca']);
                plot3(temp(1),temp(2),temp(3), ...
                    '.','MarkerSize',12,'Color',plot_colors(targ,:));
            end
            axis square;
            set(gca,'Box','off','TickDir','out','FontSize',14);
            %         xlabel('Mode 1');
            %         ylabel('Mode 2');
            %         zlabel('Mode 3');
            title([array '-' epoch_names{j} '-' td_names{i}]);
            switch lower(array)
                case 'm1'
                    set(gca,'View', [4.0000  -80.4000])
                case 'pmd'
                    set(gca,'View', [6.8000   -3.6000]);
            end
        end
    end
    
    
    linkprop(ax,'View');
    
    if save_figs
        fn = [monkey '_COvsRT_' array '_Clouds'];
        saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot kinematics
f = figure;
hold all;
bar(corr_kin.r);

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',1:2,'XTickLabel',{'X Vel','Y Vel'});
title('Kinematics correlation - CO vs RT')
set(gca,'YLim',[0 1],'XLim',[0 3])

if save_figs
    fn = [monkey '_COvsRT_BehaviorCorrelation'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot alignment results
sig_prctile = [2.5 97.5];
for iArray = 1:length(arrays)
    array = arrays{iArray};
    
    results = array_results.(array);
    
    f = figure;
    hold all;
    m = mean(results);
    s =  std(results);
    [hbar,herr] = barwitherr(s,m);
    set(hbar,'FaceColor','k','FaceAlpha',0.1,'LineWidth',2)
    set(herr,'LineWidth',2)
    
    % do statistics
    s = prctile(results(:,3)-results(:,1),sig_prctile);
    is_sig = isempty(range_intersection(s,[0 0]));
    if is_sig
        hold all;
        plot([1 3],0.98*[1 1],'k','LineWidth',2);
        text(2,0.98,'*','FontSize',30);
    end
    s = prctile(results(:,2)-results(:,1),sig_prctile);
    is_sig = isempty(range_intersection(s,[0 0]));
    if is_sig
        hold all;
        plot([1 2],0.94*[1 1],'k','LineWidth',2);
        text(1.5,0.94,'*','FontSize',30);
    end
    s = prctile(results(:,3)-results(:,2),sig_prctile);
    is_sig = isempty(range_intersection(s,[0 0]));
    if is_sig
        hold all;
        plot([2 3],0.92*[1 1],'k','LineWidth',2);
        text(2.5,0.92,'*','FontSize',30);
    end
    
    set(gca,'Box','off','TickDir','out','FontSize',14);
    set(gca,'XTick',1:3,'XTickLabel',{'Within CO','Within RT','CO to RT'});
    set(gca,'YLim',[0 1]);
    
    ylabel('Canon. Corr.');
    title(array);
    
    if save_figs
        fn = [monkey '_COvsRT_' array 'CCs'];
        saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
        saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
    end
end

% now do the normalized version
m1 = array_results.M1;
pmd = array_results.PMd;

m1 =  m1./mean(m1(:,1));
pmd =  pmd./mean(pmd(:,1));
m1 = m1(:,2:end);
pmd = pmd(:,2:end);
%
f = figure;
hold all;
m = [mean(m1); mean(pmd)];
s =  [std(m1); std(pmd)];
[hbar,herr] = barwitherr(s',m');
set(hbar,'FaceAlpha',0.3,'LineWidth',2)
set(herr,'LineWidth',2)

% do stats
s = prctile(m1(:,1)-pmd(:,1),sig_prctile);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([0.9 1.1],0.92*[1 1],'k','LineWidth',2);
    text(1,0.92,'*','FontSize',30);
end
s = prctile(m1(:,2)-pmd(:,2),sig_prctile);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([1.9 2.1],0.92*[1 1],'k','LineWidth',2);
    text(2,0.92,'*','FontSize',30);
end
s = prctile(m1(:,1)-m1(:,2),sig_prctile);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([0.9 1.9],0.96*[1 1],'k','LineWidth',2);
    text(1.4,0.96,'*','FontSize',30);
end
s = prctile(pmd(:,1)-pmd(:,2),sig_prctile);
is_sig = isempty(range_intersection(s,[0 0]));
if is_sig
    hold all;
    plot([1.1 2.1],0.96*[1 1],'k','LineWidth',2);
    text(1.6,0.96,'*','FontSize',30);
end


set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',1:3,'XTickLabel',{'Within RT','CO to RT'});
set(gca,'YLim',[0 1]);

h = legend({'M1','PMd'},'Location','South');
set(h,'Box','off','FontSize',14);

ylabel('Canon. Corr. (Norm. by Within CO)');
title(array);

if save_figs
    fn = [monkey '_COvsRT_CCsNorm_Bar'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end


% now do the histogram version
bins = 0:0.005:1;
f = figure;
hold all;
[n,x] = hist(m1(:,1),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','b','EdgeColor','none');
[n,x] = hist(pmd(:,1),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','r','EdgeColor','none');

[n,x] = hist(m1(:,2),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','c','EdgeColor','none');
[n,x] = hist(pmd(:,2),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','m','EdgeColor','none');

set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Normalized CC');
ylabel('Density');

h = legend({'M1 Within','PMd Within','M1 CO to RT','PMd CO to RT'},'Location','NorthWest');
set(h,'Box','off');

if save_figs
    fn = [monkey '_COvsRT_CCsNorm_HistWithAll'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end


f = figure('Position',[100 100 800 400]);
subplot(1,2,1);
hold all;
[n,x] = hist(m1(:,1),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','b','EdgeColor','none');
[n,x] = hist(m1(:,2),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','c','EdgeColor','none');
h = legend({'M1 Within','M1 CO to RT'},'Location','NorthWest');
set(h,'Box','off');

set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Normalized CC');
ylabel('Density');
set(gca,'XLim',[0.6 1]);


subplot(1,2,2);
hold all;
[n,x] = hist(pmd(:,1),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','r','EdgeColor','none');
[n,x] = hist(pmd(:,2),bins);
n = n./sum(n);
h = bar(x,n,'hist');
set(h,'FaceAlpha',0.3,'FaceColor','m','EdgeColor','none');

h = legend({'PMd Within','PMd CO to RT'},'Location','NorthWest');
set(h,'Box','off');

set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Normalized CC');
ylabel('Density');

set(gca,'XLim',[0.6 1]);

if save_figs
    fn = [monkey '_COvsRT_CCsNorm_HistByArray'];
    saveas(f,fullfile(save_dir,folder_name,[fn '.fig']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.pdf']));
    saveas(f,fullfile(save_dir,folder_name,[fn '.png']));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
