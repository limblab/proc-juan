%
% I THINK THE CC PART OF THE CODE NEEDS TO BE CHANGED TO MAKE IT ONLY TAKE
% THE SUBSET OF SESSIONS WE'RE INTERESTED IN
% 


% load parameters
if ~exist('proj_params','var')
    params      = batch_compare_manifold_projs_defaults;
    disp('Loading default compare_manifold_projs params');
else
    params      = proj_params;
end

% the paper was done with n = 12 instead of 20
params.dim_manifold = 12;
dims            = 1:params.dim_manifold;


% -------------------------------------------------------------------------
% choose datasets
ds              = [4:6 10:11];% [1:3 7:9] for 'Wrist' [4:6 10:11] for 'Reach-to-grasp'

% choose the subsets of parameters
params_orig     = params; % for later
    
params.time_win = params.time_win(ds,:);


% -------------------------------------------------------------------------
% do for all the trials
target          = 'all_conc'; % Not used yet. 'all_conc' by default
% do one plot per session?
plot_p_session  = false;
% sort neurons when plotting?
sort_neurons    = false;
% take abs value corrs?
abs_corr        = false;

% -------------------------------------------------------------------------
% Do for all sessions

for s = 1:length(ds)

    
    % Retrieve time window
    time_win    = params.time_win(s,:);


    % 1) equalize trial duration across all tasks
    tda         = equalize_single_trial_dur( datasets{ds(s)}.stdata, ...
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

    
    % define variables 
    corr_n      = zeros( numel(datasets{ds(s)}.neural_chs), size(comb_tds,1) );
    cc_lv       = zeros( length(dims), size(comb_tds,1) );
    
    
    % -----------------------------------------------------------------
    % do for all pairs of tasks
    
    for p = 1:size(comb_tds,1)
    
        % -----------------------------------------------------------------
        % Canonical correlation between latent variables
        
        lv1         = tda{comb_tds(p,1)}.target{end}.neural_data.dim_red.scores(:,dims);
        lv2         = tda{comb_tds(p,2)}.target{end}.neural_data.dim_red.scores(:,dims);

        [~, ~, cc]  = canoncorr( lv1, lv2 );

        cc_lv(:,p)  = cc;


        % -----------------------------------------------------------------
        % Pairwise correlations between neurons
        
        n1          = tda{comb_tds(p,1)}.target{end}.neural_data.conc_smoothed_fr;
        n2          = tda{comb_tds(p,2)}.target{end}.neural_data.conc_smoothed_fr;
        
        if abs_corr
            cr      = abs(calc_r( n1, n2 ));
        else
            cr      = calc_r( n1, n2 );
        end

        corr_n(:,p) = cr';
    end
    
    
    % ---------------------------------------------------------------------
    % STORE RESULTS IN A GLOBAL BAR
    res{s}.cc_lv    = cc_lv;
    res{s}.corr_n   = corr_n;
    
    
    % ---------------------------------------------------------------------
    % PLOT PER SESSION
    if plot_p_session
        
        % create legend
        lgnd = cell(size(comb_tds,1),1);
        for p = 1:size(lgnd,1)
            task1   = datasets{ds(s)}.labels{comb_tds(p,1)};
            task2   = datasets{ds(s)}.labels{comb_tds(p,2)};
            lgnd{p} = [task1 ' vs ' task2];
        end
        res{s}.lgnd = lgnd;
        
        % Colormap-like figure
        figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
        subplot(121)
        imagesc(cc_lv)
        caxis([0 1]), colorbar
        set(gca,'TickDir','out','FontSize',14), box off,
        set(gca,'XTick',1:6), set(gca,'XTickLabel',lgnd), set(gca,'XTickLabelRotation',45)
        ylabel('Neural mode dynamics')
        
        subplot(122)
        if sort_neurons
            imagesc(sort(corr_n,1,'descend'))
        else
            imagesc(corr_n)
        end
        caxis([0 1]), colorbar
        set(gca,'TickDir','out','FontSize',14), box off,
        set(gca,'XTick',1:6), set(gca,'XTickLabel',lgnd), set(gca,'XTickLabelRotation',45)
        ylabel('Neural mode dynamics')
        
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
        ylim([0 1]), xlim([0 dims(end)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('Correlation neural mode dynamics'), xlabel('Manifold projection')
        
        subplot(122), hold on
        for t = 1:size(cc_lv,2)
            if sort_neurons
                plot(sort(corr_n(:,t),1,'descend'),'color',cols_cc(t,:),'linewidth',1.5)
            else
                plot(corr_n(:,t),'color',cols_cc(t,:),'linewidth',1.5)
            end
        end
        ylim([0 1]), xlim([0 size(corr_n,1)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('Correlation neural units'), xlabel('Neural unit')
    end
        


end


%% ------------------------------------------------------------------------
% SUMMARY STATISTICS


% Statistics neurons
all_corrs_n         = [];
if abs_corr
    x_ax_hist       = 0:0.025:1;
else
    x_ax_hist       = -1:0.05:1;
end


for s = 1:length(ds)
    all_corrs_n     = [all_corrs_n, reshape(res{s}.corr_n,[],1)'];
end

mn_corr_n           = mean(all_corrs_n);
sd_corr_n           = std(all_corrs_n);
sem_corr_n          = sd_corr_n/sqrt(numel(all_corrs_n));

[y_hist_corr_n, x_hist_corr_n] = histcounts(all_corrs_n,x_ax_hist);


% Statistics canonical correlations
cc_lv_mtrx          = cell2mat( cellfun( @(x) x.cc_lv, res, 'UniformOutput', false ) );

% keep only the CCs that are significant (if the number of significant CCs
% are available)
if exist('proj_results','var')
    ccs_to_keep     = cell2mat( ...
                        arrayfun( @(x) x.nbr_canon_corrs_above_chance, ...
                        proj_results.summary_data.canon_corrs, ...
                        'UniformOutput', false ) );
                    
    all_cc_lv       = [];
    for c = 1:size(cc_lv_mtrx,2) % NEED TO FIX THIS !!!!!
        all_cc_lv   = [all_cc_lv, cc_lv_mtrx(1:ccs_to_keep(c),c)'];
    end
% or keep them all
else    
    all_cc_lv       = reshape(cc_lv_mtrx,[],1)'; 
end

mn_cc_lv            = mean(all_cc_lv);
sd_cc_lv            = std(all_cc_lv);
sem_cc_lv           = sd_cc_lv/sqrt(numel(all_cc_lv));


[y_hist_cc_lv, x_hist_cc_lv] = histcounts(all_cc_lv,x_ax_hist);


% normalize the counts to percentages
norm_y_hist_corr_n  = y_hist_corr_n/numel(all_corrs_n);
norm_y_hist_cc_lv   = y_hist_cc_lv/numel(all_cc_lv);



%% ------------------------------------------------------------------------
% SUMMARY PLOTS


% PLOT HISTOGRAM CORRELATION NEURONS
% 
% figure, hold on
% tb = bar(x_hist_corr_n(1:end-1),y_hist_corr_n,'histc');
% set(tb,'FaceColor','k','EdgeColor','k')
% set(gca,'TickDir','out','FontSize',14), 
% box off,
% if abs_corr
%     xlim([-.05 1.05])
% else
%     xlim([-1.05 1.05])
% end
% ylim([0 375])
% xlabel('Correlation'),ylabel('Counts')
% 
% y_stats_corr        = 350*ones(1,2);
% plot([mn_corr_n-sd_corr_n, mn_corr_n+sd_corr_n],y_stats_corr,'k','linewidth',2)
% plot(mn_corr_n,y_stats_corr(1),'.k','markersize',28)
% 

% PLOT NORMALIZED HISTOGRAM CORRELATION NEURONS

figure, hold on
tb = bar(x_hist_corr_n(1:end-1),y_hist_corr_n/numel(all_corrs_n),'histc');
set(tb,'FaceColor','k','EdgeColor','k')
set(gca,'TickDir','out','FontSize',14), 
box off,
if abs_corr
    xlim([-.05 1.05])
else
    xlim([-1.05 1.05])
end
ylim([0 .2])
xlabel('Correlation neural activity across tasks'),ylabel('Normalized Counts (%)')

plot([mn_corr_n-sd_corr_n, mn_corr_n+sd_corr_n],[.19 .19],'k','linewidth',2)
plot(mn_corr_n,.19,'.k','markersize',28)

text(-.9,.18,['n=' num2str(length(all_corrs_n))],'FontSize',14);


% PLOT HISTOGRAM CORRELATION NEURONS VS CANONICAL CORR LATENT VARS
% 
% figure, hold on
% tb = bar(x_hist_corr_n(1:end-1),y_hist_corr_n,'histc');
% set(tb,'FaceColor','k','EdgeColor','k')
% tb2 = bar(x_hist_cc_lv(1:end-1),y_hist_cc_lv,'histc');
% tb2.FaceAlpha = 0.5;
% set(tb2,'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
% set(gca,'TickDir','out','FontSize',14), 
% box off,
% legend('Neurons','Neural modes','Location','NorthWest'), legend boxoff
% if abs_corr
%     xlim([-.05 1.05])
% else
%     xlim([-1.05 1.05])
% end
% ylim([0 400])
% xlabel('Correlation'),ylabel('Counts')
% 
% % Add some stats to the plot
% y_stats_corr        = 350*ones(1,2);
% y_stats_cc          = 375*ones(1,2);
% plot([mn_corr_n-sd_corr_n, mn_corr_n+sd_corr_n],y_stats_corr,'k','linewidth',2)
% plot(mn_corr_n,y_stats_corr(1),'.k','markersize',28)
% plot([mn_cc_lv-sd_cc_lv, mn_cc_lv+sd_cc_lv],y_stats_cc,'color',[.6 .6 .6],'linewidth',2)
% plot(mn_cc_lv,y_stats_cc(1),'.','color',[.6 .6 .6],'markersize',28)
% ylim([0 400])


% % PLOT NORMALIZED VERSION PREVIOUS HISTOGRAM
% 
% figure, hold on
% tb = bar(x_hist_corr_n(1:end-1),norm_y_hist_corr_n,'histc');
% set(tb,'FaceColor','k','EdgeColor','k')
% tb2 = bar(x_hist_cc_lv(1:end-1),norm_y_hist_cc_lv,'histc');
% set(tb2,'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
% tb2.FaceAlpha = 0.5; tb2.EdgeAlpha = 0.5;
% set(gca,'TickDir','out','FontSize',14)
% box off,
% legend('Neurons','Neural modes','Location','NorthWest'), legend boxoff
% if abs_corr
%     xlim([-.05 1.05])
% else
%     xlim([-1.05 1.05])
% end
% ylim([0 0.18])
% xlabel('Correlation'),ylabel('Normalized counts')
% 
% % Add some stats to the plot
% y_stats_corr        = 0.18*ones(1,2);
% y_stats_cc          = 0.19*ones(1,2);
% plot([mn_corr_n-sd_corr_n, mn_corr_n+sd_corr_n],y_stats_corr,'k','linewidth',2)
% plot(mn_corr_n,y_stats_corr(1),'.k','markersize',28)
% plot([mn_cc_lv-sd_cc_lv, mn_cc_lv+sd_cc_lv],y_stats_cc,'color',[.6 .6 .6],'linewidth',2)
% plot(mn_cc_lv,y_stats_cc(1),'.','color',[.6 .6 .6],'markersize',28)
% ylim([0 0.20])



% save all correlations
if ~exist('all_r','var')
    all_r{1} = all_corrs_n;
else
    all_r{length(all_r)+1} = all_corrs_n;
end

% Clear vars
params = params_orig;
clearvars -except *_results *_params datasets manifold_dim params all_r