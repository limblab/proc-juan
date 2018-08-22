%
%
% Comparison of manifold orientation dynamics using single units and
% threshold crossings
%


plot_ISIs_yn                = true;

pathBDF = '/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/sorted unit dataset';

cd(pathBDF)
mani_dims                   = 1:12;
mu_ds                       = 8; % multiunit dataset

% Pre-process the single unit datasets BDFs
prep_params                 = batch_preprocess_dim_red_data_params_defaults;
prep_params.bin_size        = .02;
prep_params.save_data       = false;
prep_params.norm_trial_data = 'min_dur'; % 'min_dur' 't_warp'
prep_params.dim_red_emg     = 'nmf'; % 'nmf'
prep_params.emg_factors     = 3;

proj_params                 = batch_compare_manifold_projs_defaults;


% preprocess the data
sortdatasets    = batch_preprocess_dim_red_data( pathBDF, prep_params );
sortdatasets    = sortdatasets{1};

% load all the data
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% just keep the same session that we have sorted
mudatasets      = datasets{mu_ds};


% COMPARE STUFF !!!
% Do CCA between the latent activities based on sorted units and
% multi-units
n_ds            = length(sortdatasets.labels);
ccs             = zeros(length(mani_dims),n_ds);
pas             = zeros(length(mani_dims),n_ds);

% equalize single trial duration 
mustdata        = mudatasets.stdata;
sortstdata      = sortdatasets.stdata;

% 1) equalize trial duration across all tasks
mustdata = equalize_single_trial_dur( mustdata, 'time_win', proj_params.time_win(mu_ds,:) );
sortstdata = equalize_single_trial_dur( sortstdata, 'time_win', proj_params.time_win(mu_ds,:) );
    
% 2) equalize number of trials for all targets of a given task
for i = 1:length(mustdata)
    mustdata{i} = equalize_nbr_trials_p_target( mustdata{i} );
    sortstdata{i} = equalize_nbr_trials_p_target( sortstdata{i} );
end
    
% 3) equalize number of trials across tasks
mustdata = equalize_nbr_trials_across_tasks( mustdata, 'all_conc' );
sortstdata = equalize_nbr_trials_across_tasks( sortstdata, 'all_conc' );


for d = 1:n_ds
    
    % CCA -- we need
    lv_mu   = mustdata{d}.target{end}.neural_data.dim_red.scores(:,mani_dims);
    lv_sort = sortstdata{d}.target{end}.neural_data.dim_red.scores(:,mani_dims);
    
    [~,~,ccs(:,d)] = canoncorr(lv_mu,lv_sort);
    
%     % PAs
%     w_mu        = mudatasets.dim_red_FR{d}.w(:,mani_dims);
%     w_sort      = sortdatasets.dim_red_FR{d}.w(:,mani_dims);
%     
%     % If the number of multiunits and sorted units is different take a subset
%     % of the multiunits
%     if size(w_mu,1) ~= size(w_sort,1)
%         if size(w_mu,1) > size(w_sort,1)
%             
%             fr_mus = mustdata{d}.target{end}.neural_data.conc_smoothed_fr;
%             fr_mus = fr_mus(:,randperm(size(w_mu,1)));
%             fr_mus = fr_mus(:,1:size(w_sort,1));
%             
%             w_mu = pca(fr_mus);
%         else
%            error('not implemented yet'); 
%         end
%     end
%     
%     pas(:,d)    = principal_angles(w_mu,w_sort);
end


% PLOT COMPARISONS manifold and latent activity comparisons
figure,
% subplot(121)
% plot(rad2deg(pas),'linewidth',1.5)
% ylim([0 1]),xlim([0 mani_dims(end)])
% set(gca,'FontSize',14,'TickDir','out'); box off
% legend(mudatasets.labels),legend boxoff
% xlabel('Neural mode'); ylabel('Principal angle (deg)')
% subplot(122)
plot(ccs,'linewidth',1.5)
ylim([0 1]),xlim([0 mani_dims(end)])
set(gca,'FontSize',14,'TickDir','out'); box off
set(gcf,'color',[1 1 1])
xlabel('Neural mode'); ylabel('CC latent activity')


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Comparison of manifolds across tasks

comb_tasks = nchoosek(1:4,2);
n_comb_tasks = size(comb_tasks,1);


cross_task_pas = zeros(length(mani_dims),n_comb_tasks);
cross_task_ccs = zeros(length(mani_dims),n_comb_tasks);

for c = 1:n_comb_tasks
   
    t1 = comb_tasks(c,1);
    t2 = comb_tasks(c,2);

    % For PRINCIPAL ANGLES
    w1 = sortdatasets.dim_red_FR{t1}.w(:,mani_dims);
    w2 = sortdatasets.dim_red_FR{t2}.w(:,mani_dims);
    
%     % For CCs
%     lv1 = sortdatasets.dim_red_FR{t1}.scores(:,mani_dims);
%     lv2 = sortdatasets.dim_red_FR{t2}.scores(:,mani_dims);
    
    % check that the dimensionality of the space is the same for both tasks 
    % -- in one of them the cov matrix is not full rank  
    if size(w1,1) ~= size(w2,1)
        
        fr1 = sortstdata{t1}.target{end}.neural_data.conc_smoothed_fr;
        fr2 = sortstdata{t2}.target{end}.neural_data.conc_smoothed_fr;

        if size(w1,1) < size(w2,1)

            % diff_dims = size(w2,1) - size(w1,1);
            rnd_units = randperm(size(fr1,2));

        elseif size(w2,1) < size(w1,1)

            % diff_dims = size(w1,1) - size(w2,1);
            rnd_units = randperm(size(fr2,2));
        end
        
        units_keep = rnd_units;
        
        fr1 = fr1(:,units_keep);
        fr2 = fr2(:,units_keep);
        
        [w1, lv1] = pca(fr1);
        w1 = w1(:,mani_dims);
        lv1 = lv1(:,mani_dims);
        
        [w2, lv2] = pca(fr2);
        w2 = w2(:,mani_dims);
        lv2 = lv2(:,mani_dims);
    end
    
    cross_task_pas(:,c) = principal_angles(w1,w2);
    
%     [~,~,cross_task_ccs(:,c)] = canoncorr(lv1,lv2);
end

% load significance threshold
if ~exist('angle_non_orth','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets_original_submission.mat');
end

signif_th = angle_non_orth(:,1,8);

lgnd = {};
for c = 1:n_comb_tasks
    t1 = comb_tasks(c,1);
    t2 = comb_tasks(c,2);
    lgnd{c} = [sortdatasets.labels{t1} ' vs ' sortdatasets.labels{t2}];
end

figure, 
hold on
plot(rad2deg(cross_task_pas),'linewidth',1.5)
plot(signif_th,'linewidth',1.5,'linestyle','-.','color',[.65 .65 .65])
set(gca,'FontSize',14,'TickDir','out'); box off
% legend(mudatasets.labels),legend boxoff
ylim([0 90]),xlim([0 mani_dims(end)])
legend(lgnd,'location','SouthEast'), legend boxoff
xlabel('Neural mode'); ylabel('Principal angle (deg)')
set(gcf,'color','w')


% subplot(122)
% plot(rad2deg(cross_task_ccs),'linewidth',1.5)
% ylim([0 1]),xlim([0 mani_dims(end)])
% set(gca,'FontSize',14,'TickDir','out'); box off
% set(gcf,'color',[1 1 1])
% xlabel('Neural mode'); ylabel('CC latent activity')


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Plot a few ISIs and their mean waveforms

if plot_ISIs_yn

    if ~exist('cbdf','var')
        load([pathBDF filesep 'all_tasks_Jango_20140923_sorted.mat'])
    end
    
    units = [1:4 10 11 13 15:17 18 19];

    n_rows = floor(sqrt(length(units)));
    n_cols = ceil(length(units)/n_rows);

    x_hist = 0:0.005:1.005;

    cols   = parula(length(sortstdata)+1);

    
    % ISIs
    figure
    for u = 1:length(units)

        subplot(n_rows,n_cols,u), hold on

        % for each task
        for t = 1:length(sortstdata)

            ucomp = sortdatasets.binned_data{t}.neuronIDs == neural_chs(units(u),:);
            idx_unit = find(sum(ucomp,2)==2);

            spiking = cbdf(t).units(idx_unit).ts;
            t_y = histcounts(diff(spiking),x_hist)/numel(spiking);
            plot(x_hist(1:end-1),t_y,'color',cols(t,:),'linewidth',1.5);
        end

        xlim([0 0.300])
        title(num2str(neural_chs(idx_unit,:)))
        set(gca,'FontSize',14,'TickDir','out'); box off
    end
    subplot(n_rows,n_cols,1),legend(sortdatasets.labels), legend boxoff
    set(gcf,'color','w')
    
    
    % MEAN WAVEFORMS
    figure
    for u = 1:length(units)
       
        subplot(n_rows,n_cols,u), hold on
        
        % for each task
        for t = 1:length(sortstdata)
            
            ucomp = sortdatasets.binned_data{t}.neuronIDs == neural_chs(units(u),:);
            idx_unit = find(sum(ucomp,2)==2);
            
            
            waveforms = cbdf(t).units(idx_unit).waveforms;
            mn_w = mean(waveforms,1);
            t_w = 0:1/30000:1/30000*(length(mn_w)-1);
            plot(t_w,mn_w,'color',cols(t,:),'linewidth',1.5);
        end
        
        xlim([0 1/30000*(length(mn_w)-1)])
        title(num2str(neural_chs(idx_unit,:)))
        set(gca,'FontSize',14,'TickDir','out'); box off
    end
    subplot(n_rows,n_cols,1),legend(sortdatasets.labels,'Location','SouthEast'), legend boxoff
    set(gcf,'color','w')
end
