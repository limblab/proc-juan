%
% Canonical correlation between EMGs
%



% load parameters
params          = batch_compare_manifold_projs_defaults;
% do for all the trials
target          = 'all_conc';
% do one plot per session?
plot_p_session  = false;
% sort neurons when plotting?
sort_emgs       = false;
% take abs value corrs?
abs_corr        = false;
% compare to CC neural mode dynamics
comp_cc_lvs     = true;

plot_pairs_diff_cols = false;


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

    % get all pairwise combinations of days
    comb_tds    = nchoosek(1:length(tda),2);

    
    % define variables 

    % only use the "chosen EMGs"
    nbr_ch_emgs = size(tda{comb_tds(1,1)}.target{end}.emg_data.conc_emg,2);
    % preallocate matrices
    cc_emg      = zeros( nbr_ch_emgs, size(comb_tds,1) );
    r_emg       = zeros( nbr_ch_emgs, size(comb_tds,1) );
    
    
    % -----------------------------------------------------------------
    % do for all pairs of tasks
    
    for p = 1:size(comb_tds,1)
    
        % -----------------------------------------------------------------
        % Canonical correlation between latent variables
        
        emg1        = tda{comb_tds(p,1)}.target{end}.emg_data.conc_emg;
        emg2        = tda{comb_tds(p,2)}.target{end}.emg_data.conc_emg;

        [~, ~, cc]  = canoncorr( emg1, emg2 );

        cc_emg(:,p)  = cc;

        % -----------------------------------------------------------------
        % Pairwise correlations between EMGs, to assess significance
        
        if abs_corr
            cr      = abs(calc_r( emg1, emg2 ));
        else
            cr      = calc_r( emg1, emg2 );
        end

        r_emg(:,p) = cr';
    end
    
    
    % ---------------------------------------------------------------------
    % STORE RESULTS IN A GLOBAL BAR
    res{s}.cc_emg   = cc_emg;
    res{s}.r_emg    = r_emg;
   
    
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
        
        % Plots with traces
        if size(cc_emg,2) > 1
            cols_cc = parula(size(cc_emg,2));
        else
            cols_cc = [0 0 0];
        end
        figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
        subplot(121), hold on
        for t = 1:size(cc_emg,2)
            plot(cc_emg(:,t),'color',cols_cc(t,:),'linewidth',1.5)
        end
        ylim([0 1]), xlim([0 nbr_ch_emgs(end)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('Correlation neural mode dynamics'), xlabel('Manifold projection')
        
        subplot(122), hold on
        for t = 1:size(cc_emg,2)
            if sort_emgs
                plot(sort(r_emg(:,t),1,'descend'),'color',cols_cc(t,:),'linewidth',1.5)
            else
                plot(r_emg(:,t),'color',cols_cc(t,:),'linewidth',1.5)
            end
        end
        ylim([0 1]), xlim([0 nbr_ch_emgs(end)])
        set(gca,'TickDir','out','FontSize',14), box off,
        legend(lgnd), legend boxoff
        ylabel('Correlation neural units'), xlabel('Neural unit')
    end
end

%% ------------------------------------------------------------------------
%   SUMMARY STATS?

all_cc_emg          = [];

for s = 1:length(datasets)
    all_cc_emg      = [all_cc_emg, reshape(res{s}.cc_emg,1,[])];
end

mn_cc_emg           = mean(all_cc_emg);
sd_cc_emg           = std(all_cc_emg);
sem_cc_emg          = sd_cc_emg/sqrt(numel(all_cc_emg));



%% ------------------------------------------------------------------------
% COMPARE CC NEURAL MODE DYNAMICS

if comp_cc_lvs
   
    % load paper version of all the results
    load('/Users/juangallego/Documents/Publications/2017 - Multi task manifold + some stuff not on DB/Raw data for figs - 2017-06-20/Results_manifold_notebook_v1_20170620.mat');
    
    % get metadata
    metad           = batch_get_monkey_task_data(datasets);
    diff_cols       = parula(length(metad.task_pairs.unique_pairs));
    % ctr for colors
    ctr_cols        = 0;
    
    % variable to store the neural mode CC vs EMG CC
    cc_ratio_p      = [];
    
    
    % --------------------------------------------------------------------    
    % plot the neural CC to EMG CC ratio, for each task comparison
%     figure, hold on
    for s = 1:length(datasets)
    
       cc_lv        = proj_results.data{s}.can_corrs.cc';
       cc_emg       = res{s}.cc_emg;
       
       % keep only as many neural CCs as EMG CCs
       only_cc_lv   = cc_lv(1:size(cc_emg,1),:);
       
       % compute the ratio of the neural CCs to the EMG CCs
       cc_ratio     = only_cc_lv./cc_emg;
       
       % pool over sessions
       cc_ratio_p   = [cc_ratio_p, reshape(cc_ratio,1,[])];
       
%        % figure       
%        if ~plot_pairs_diff_cols
%            plot(cc_ratio,'k','linewidth',1)
%        else
%            % increase ctr for colors
%            for p = 1:size(only_cc_lv,2)
%                
%                % COMPROBAR !!!!
%                tc = ctr_cols + p;
%                plot(cc_ratio(:,p),'k','linewidth',1,'color', diff_cols(metad.task_pairs.task_pair_nbr(tc),:),'linewidth',1.5)
%                % COMPROBAR !!!!
%                
%            end
%            ctr_cols = ctr_cols + size(only_cc_lv,2);
%        end
    end
%     plot([1 length(proj_results.data{s}.can_corrs.cc)],[1 1],...
%             'linewidth',6,'color',[.6 .6 .6],'linestyle','--');
%     ylabel('Neural mode dynamics vs EMG correlation')
%     xlabel('Projection'), ylim([0 12])
%     set(gca,'TickDir','out','FontSize',14), box off,
%     xlabel('Projection'),ylabel('Neural mode dynamics to EMG Canon. Corr. ratio')

    
    % --------------------------------------------------------------------
    % plot all the Neural CCs and EMG CCs for each session
    if plot_p_session
        for s = 1:length(datasets)
            figure,hold on
            plot(proj_results.data{s}.can_corrs.cc','k','linewidth',1.5)
            plot(res{s}.cc_emg,'color',[.6 .6 .6],'linewidth',1.5)
            ylim([0 1]),xlim([0 length(proj_results.data{s}.can_corrs.cc)])
            set(gca,'TickDir','out','FontSize',14), box off,
            xlabel('Projection'),ylabel('Canonical correlation')
        end
    end
    % --------------------------------------------------------------------
    % histogram of the Neural CC to EMG CC ratio

    % remove outliers
    cc_ratio_p_noout = cc_ratio_p(cc_ratio_p<prctile(cc_ratio_p,95));
    
    mn_cc_ratio     = mean(cc_ratio_p_noout);
    sd_cc_ratio     = std(cc_ratio_p_noout);
    
    % compute histogram
    x_hist_cc_ratio = 0:0.5:12;
    [y_hist_cc_ratio, x_hist_cc_ratio] = histcounts( cc_ratio_p, x_hist_cc_ratio );
    
    y_stats         = max(y_hist_cc_ratio);
    
%     % plot
%     figure, hold on
%     tb = bar(x_hist_cc_ratio(1:end-1),y_hist_cc_ratio,'histc');
%     set(tb,'FaceColor','k','EdgeColor','k')
%     set(gca,'TickDir','out','FontSize',14), box off,
%     plot([1 1],[0 max(y_hist_cc_ratio)],'linewidth',6,'color',[.6 .6 .6],'linestyle','--')
%     plot([mn_cc_ratio-sd_cc_ratio,mn_cc_ratio+sd_cc_ratio],[y_stats+2.5 y_stats+2.5],...
%         'k','linewidth',2);
%     plot(mn_cc_ratio,y_stats+2.5,'.k','markersize',28);
%     ylabel('Counts'),xlabel('Neural mode dynamics to EMG Canon. Corr. ratio')
%     ylim([0 max(y_hist_cc_ratio)+5])
    
    % plot normalized 
    figure, hold on
    tb = bar(x_hist_cc_ratio(1:end-1),y_hist_cc_ratio/numel(cc_ratio_p),'histc');
    set(tb,'FaceColor','k','EdgeColor','k')
    set(gca,'TickDir','out','FontSize',14), box off,
    plot([1 1],[0 max(y_hist_cc_ratio)],'linewidth',6,'color',[.6 .6 .6],'linestyle','--')
    plot([mn_cc_ratio-sd_cc_ratio,mn_cc_ratio+sd_cc_ratio],[y_stats y_stats]/numel(cc_ratio_p)+.01,...
        'k','linewidth',2);
    plot(mn_cc_ratio,y_stats/numel(cc_ratio_p)+.01,'.k','markersize',28);
    ylabel('Normalized Counts'),xlabel('Neural mode dynamics to EMG Canon. Corr. ratio')
    ylim([0 max(y_hist_cc_ratio/numel(cc_ratio_p))+.02])
    text(10,0.2,['n=' num2str(numel(cc_ratio_p))],'FontSize',14)
end



