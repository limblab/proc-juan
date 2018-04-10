

% load parameters
if ~exist('proj_params','var')
    params      = batch_compare_manifold_projs_defaults;
    disp('Loading default compare_manifold_projs params');
else
    params      = proj_params;
end

% the paper was done with n = 12 instead of 20
params.dim_manifold = 12;
dims            = 'same_as_emg'; % 'same_as_emg'

% do for all the trials
target          = 'all_conc'; % Not used yet. 'all_conc' by default
% do one plot per session?
plot_p_session  = false;



% -------------------------------------------------------------------------
% Do for all sessions

for s = 1:length(datasets)

    
    % Retrieve time window
    time_win    = params.time_win(s,:);


    % 1) equalize trial duration across all tasks
    tda         = equalize_single_trial_dur( datasets{s}.stdata, ...
                    'time_win', time_win );

    % 2) equalize number of trials for all targets of a given task
    for i = 1:numel(tda)
        tda{i}  = equalize_nbr_trials_p_target( tda{i} );
    end

    % 3) equalize number of trials across tasks
    tda         = equalize_nbr_trials_across_tasks( tda, target );



    % ---------------------------------------------------------------------
    % do the analysis


    clear canon_corr;

    % get all pairwise combinations of tasks
    comb_tds    = nchoosek(1:length(tda),2);

    % Nbr of EMGs
    n_emgs      = length(tda{1}.target{1}.emg_data.emg_names);
    
    % See if we want the neural manifold dimensionality == EMG manifold
    % dimensionality 
    if ischar(dims)
       if strcmp(dims,'same_as_emg')
           t_dim = 1:n_emgs;
       else
          error('DIMS has to be a scalar or same_as_emg');
       end
    else
        t_dim = dims;
    end
    
    cc_lv           = zeros( length(t_dim), size(comb_tds,1) );
    cc_emgs         = zeros( n_emgs, size(comb_tds,1) );
    
    
    % -----------------------------------------------------------------
    % do for all pairs of tasks
    
    for p = 1:size(comb_tds,1)
    
        % -----------------------------------------------------------------
        % Canonical correlation between latent variables
        
        lv1         = tda{comb_tds(p,1)}.target{end}.neural_data.dim_red.scores(:,t_dim);
        lv2         = tda{comb_tds(p,2)}.target{end}.neural_data.dim_red.scores(:,t_dim);

        [~, ~, cc]  = canoncorr( lv1, lv2 );

        cc_lv(:,p)  = cc;


        % -----------------------------------------------------------------
        % Canonical correlation between EMGs
        
        emg1        = tda{comb_tds(p,1)}.target{end}.emg_data.conc_emg;
        emg2        = tda{comb_tds(p,2)}.target{end}.emg_data.conc_emg;
        
        [~, ~, cc_e] = canoncorr( emg1, emg2 );
        
        cc_emgs(:,p) = cc_e;
    end
    
    
    % ---------------------------------------------------------------------
    % STORE RESULTS IN A GLOBAL BAR
    res{s}.cc_lv    = cc_lv;
    res{s}.cc_emgs  = cc_emgs;
    
    
    % ---------------------------------------------------------------------
    % PLOT PER SESSION
    if plot_p_session
        
        % create legend
        lgnd = cell(size(comb_tds,1),1);
        for p = 1:size(lgnd,1)
            task1   = datasets{s}.labels{comb_tds(p,1)};
            task2   = datasets{s}.labels{comb_tds(p,2)};
            lgnd{p} = [task1 ' vs ' task2];
        end
        res{s}.lgnd = lgnd;
        
%         % Colormap-like figure
%         figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
%         subplot(121)
%         imagesc(cc_lv)
%         caxis([0 1]), colorbar
%         set(gca,'TickDir','out','FontSize',14), box off,
%         set(gca,'XTick',1:6), set(gca,'XTickLabel',lgnd), set(gca,'XTickLabelRotation',45)
%         ylabel('Neural mode dynamics')
%         subplot(122)
%         imagesc(cc_emgs)
%         caxis([0 1]), colorbar
%         set(gca,'TickDir','out','FontSize',14), box off,
%         set(gca,'XTick',1:6), set(gca,'XTickLabel',lgnd), set(gca,'XTickLabelRotation',45)
%         ylabel('EMGs')
        
        % Plots with traces
        if size(cc_lv,2) > 1
            cols_cc = parula(size(cc_lv,2));
        else
            cols_cc = [0 0 0];
        end
        figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
        subplot(121), hold on
        for t = 1:size(cc_lv,2)
            plot(cc_lv(:,t),'color',cols_cc(t,:),'linewidth',1.5)
        end
        ylim([0 1]), xlim([0 t_dim(end)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('CC neural mode dynamics'), xlabel('Manifold projection')
        
        subplot(122), hold on
        for t = 1:size(cc_lv,2)
            plot(cc_emgs(:,t),'color',cols_cc(t,:),'linewidth',1.5)
        end
        ylim([0 1]), xlim([0 t_dim(end)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('CC EMGs'), xlabel('Neural unit')
    end
end




%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% SUMMARY STATISTICS


x_ax_hist           = 0:0.025:1;

% total number of comparisons
n_comps             = sum(cellfun(@(x) size(x.cc_lv,2), res));

% nbr of EMGs per session
n_emgs_p_ds         = cellfun(@(x) size(x.cc_emgs,1), res);

% preallocate matrices
cc_emg_mtrx         = nan(max(n_emgs_p_ds),n_comps);
if ischar(dims)
    cc_lv_mtrx      = nan(max(n_emgs_p_ds),n_comps);
else
    cc_lv_mtrx      = nan(dims(end),n_comps);
end


% -----------------------------------
% Populate data matrices

% For the neural modes
if ~ischar(dims)
    cc_lv_mtrx      = cell2mat( cellfun( @(x) x.cc_lv, res, 'UniformOutput', false ) );
else
    ptr = 1;
    for d = 1:length(datasets)
        n_comps_t   = size(res{d}.cc_lv,2);
        idxs_t      = ptr:(ptr+n_comps_t-1);
        dims_t      = 1:size(res{d}.cc_lv,1);
        cc_lv_mtrx(dims_t,idxs_t) = res{d}.cc_lv;
        ptr         = ptr + n_comps_t;
    end
end

% For the EMGs
ptr = 1;
for d = 1:length(datasets)
    n_emgs_t    = size(res{d}.cc_emgs,2);
    idxs_t      = ptr:(ptr+n_emgs_t-1);
    dims_t      = 1:size(res{d}.cc_emgs,1);
    cc_emg_mtrx(dims_t,idxs_t) = res{d}.cc_emgs;
    ptr         = ptr + n_emgs_t;
end


% -------------------------------------------------------------------------
% Summary statistics 

% % keep only the CCs that are significant (if the number of significant CCs
% % are available) -- LEGACY AND NOT IMPLEMENTED FOR THE EMGS
% if exist('proj_results','var')
%     ccs_to_keep     = cell2mat( ...
%                         arrayfun( @(x) x.nbr_canon_corrs_above_chance, ...
%                         proj_results.summary_data.canon_corrs, ...
%                         'UniformOutput', false ) );
%                     
%     all_cc_lv       = [];
%     for c = 1:size(cc_lv_mtrx,2)
%         all_cc_lv   = [all_cc_lv, cc_lv_mtrx(1:ccs_to_keep(c),c)'];
%     end
% % or keep them all
% else    
%     all_cc_lv       = reshape(cc_lv_mtrx,[],1)'; 
% end


% Reshape the data matrices to vectors without NaNs
all_cc_lv       = reshape(cc_lv_mtrx,[],1)'; 
all_cc_lv(isnan(all_cc_lv)) = [];

all_cc_emgs     = reshape(cc_emg_mtrx,[],1)'; 
all_cc_emgs(isnan(all_cc_emgs)) = [];


% For the neural modes
mn_cc_lv        = mean(all_cc_lv);
sd_cc_lv        = std(all_cc_lv);
sem_cc_lv       = sd_cc_lv/sqrt(numel(all_cc_lv));

[y_hist_cc_lv, x_hist_cc_lv] = histcounts(all_cc_lv,x_ax_hist);


% For the EMGs
mn_cc_emgs      = mean(all_cc_emgs);
sd_cc_emgs      = std(all_cc_emgs);
sem_cc_emgs     = sd_cc_emgs/sqrt(numel(all_cc_emgs));

[y_hist_cc_emgs, x_hist_cc_emgs] = histcounts(all_cc_emgs,x_ax_hist);


% normalize the counts to percentages
norm_y_hist_cc_lv   = y_hist_cc_lv/numel(all_cc_lv)*100;
norm_y_hist_cc_emgs = y_hist_cc_emgs/numel(all_cc_emgs)*100;


% -------------------------------------------------------------------------
% Ratio between EMG modes and neural modes
if length(all_cc_emgs) == length(all_cc_lv)

    hist_bins = 0.5;
    lv_to_emg_ratio = all_cc_lv./all_cc_emgs;
    x_ratio_hist = 0:hist_bins:8;
    [y_hist_ratio, x_hist_ratio] = histcounts(lv_to_emg_ratio,x_ratio_hist);
    
    norm_y_hist_ratio = y_hist_ratio/numel(all_cc_emgs)*100;
end



%% ------------------------------------------------------------------------
% SUMMARY PLOTS


% -------------------------------------------------------------------------
% PLOT NORMALIZED HISTOGRAM CC MODES AND EMGs

figure, hold on
tb = bar(x_hist_cc_emgs(1:end-1),norm_y_hist_cc_emgs,'histc');
set(tb,'FaceColor','k','EdgeColor','k')
tb2 = bar(x_hist_cc_lv(1:end-1),norm_y_hist_cc_lv,'histc');
set(tb2,'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
tb2.FaceAlpha = 0.5; tb2.EdgeAlpha = 0.5;
set(gca,'TickDir','out','FontSize',14)
box off,xlim([-.05 1.05])
xlabel('Correlation'),ylabel('Percentage (%)')
yl = ylim; ylim([0 yl(2)+3]);

% Add some stats to the plot
y_stats_emgs        = (yl(2)+1.5)*ones(1,2);
y_stats_cc          = (yl(2)+2)*ones(1,2);
plot([mn_cc_emgs-sd_cc_emgs, mn_cc_emgs+sd_cc_emgs],y_stats_emgs,'k','linewidth',2)
plot(mn_cc_emgs,y_stats_emgs(1),'.k','markersize',28)
plot([mn_cc_lv-sd_cc_lv, mn_cc_lv+sd_cc_lv],y_stats_cc,'color',[.6 .6 .6],'linewidth',2)
plot(mn_cc_lv,y_stats_cc(1),'.','color',[.6 .6 .6],'markersize',28)
ylim([0 20])
legend('EMGs','Neural modes','Location','NorthEast'), legend boxoff



% -------------------------------------------------------------------------
% PLOT RATIO CC NEURAL MODES / CC EMGs

if exist('lv_to_emg_ratio','var')
    figure,hold on
    tbr = bar(x_hist_ratio(1:end-1),norm_y_hist_ratio,'histc');
    set(tbr,'FaceColor','k','EdgeColor','k')
    yl = ylim; xlim([0 x_ratio_hist(end)]);
    plot([1 1],yl,'--','linewidth',3,'color',[.5 .5 .5])
    set(gca,'TickDir','out','FontSize',14), box off,
    xlabel('Neural mode to EMG CC ratio (same dimensionalities)')
    ylabel('Percentage (%)')
    text(x_ratio_hist(end)-2,yl(2)-5,['n=' num2str(length(lv_to_emg_ratio))],'Fontsize',14)
end

    
% -------------------------------------------------------------------------
% Clear vars
clearvars -except *_results *_params datasets manifold_dim