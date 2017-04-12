%
% COMPARE MANIFOLD PROJS OVER DAYS
%


clearvars -except bdf

%% ------------------------------------------------------------------------
% Define Parameters

dims        = 1:10;
time_win    = [0 0.7];
target      = 'all_conc';
manifold_from_day = 'same_day'; % 'same_day' 'day_1'
bdfs_to_use = [2 4:8];
demean_FRs  = true;
neural_chs  = 1:90; % Channels 67 and 70 were disabled the last two sessions
within_ccs  = false;

% relate session numbers to day between sessions
% days_btw_sessions   = [-(14+23),0,3,4,6,7,18];
days_btw_sessions   = [0 1 5 13 24:28 35];
days_btw_sessions   = days_btw_sessions(bdfs_to_use);


%% ------------------------------------------------------------------------
% Preprocess the data

% 0) Get single trial data if it doesn't exist
for i = 1:length(bdfs_to_use), 
    tda{i} = call_get_single_trial_data( bdf(bdfs_to_use(i)), 'neural_chs', neural_chs ); 
end;

% 1) equalize trial duration across all tasks
tda2        = equalize_single_trial_dur( tda, 'time_win', time_win );

% 2) equalize number of trials for all targets of a given task
for i = 1:numel(tda2)
    tda2{i} = equalize_nbr_trials_p_target( tda2{i} );
end

% 3) equalize number of trials across tasks
tda2 = equalize_nbr_trials_across_tasks( tda2, target );
    

%% ------------------------------------------------------------------------
% do the analysis

% get all pairwise combinations of days
comb_tds    = nchoosek(1:length(bdfs_to_use),2);

% predefine the matrix with plain old corrs for assessing how well CCA does
r_ctrl      = zeros(size(comb_tds,1),length(dims));


for p = 1:size(comb_tds,1)
    
    % get projected, aligned-in-time scores for this task
    switch target
        case 'all_conc'
            scores_1 = tda2{comb_tds(p,1)}.target{end}.neural_data.dim_red.scores(:,dims);
            switch manifold_from_day
                case 'day_1'
                    W       = tda2{comb_tds(p,1)}.target{end}.neural_data.dim_red.w(:,dims);
                    FRs     = tda2{comb_tds(p,2)}.target{end}.neural_data.conc_smoothed_fr;
                                        
                    % WARNING: THE SCORES WON'T BE THE SAME AS WHEN
                    % COMPUTED FROM THE ORIGINAL DATA, BECAUSE THE MEAN IS
                    % DIFFERENT (AS WE HAVE EQUALIZED THE NUMBER OF TRIALS
                    % AND TRIAL DURATION)
                    if demean_FRs
                        scores_2 = W' * ( FRs - repmat(mean(FRs)',1,size(FRs,1))' )';
                    else
                        scores_2 = W' * FRs';
                    end
                    scores_2 = scores_2';
                case 'same_day'
                    scores_2 = tda2{comb_tds(p,2)}.target{end}.neural_data.dim_red.scores(:,dims);
            end
%         otherwise
%             scores_1 = tda2{comb_tds(p,1)}.target{target}.neural_data.dim_red.scores(:,dims);
%             scores_2 = tda2{comb_tds(p,2)}.target{target}.neural_data.dim_red.scores(:,dims);            
    end
    
    
    % do CCA
   [A, B, r, U, V, stats] = canoncorr( scores_1, scores_2 ); 
   
   % do the control analysis, just correlate the scores
   for d = 1:length(dims)
       aux_r_ctrl   = corrcoef( scores_1(:,d), scores_2(:,d) );
       r_ctrl(p,d)  = aux_r_ctrl(1,2);
   end

   
   % do within task-CCs to assess significance
   if within_ccs
       cc_within_1 = canon_corr_within_task( tda{comb_tds(p,1)}, dims );
       cc_within_2 = canon_corr_within_task( tda{comb_tds(p,2)}, dims );
   
       cc_within_avg = mean([cc_within_1.cc;cc_within_2.cc]);
   end
   
   
   % store, to return
   can_corrs.lin_transform(p).A = A;
   can_corrs.lin_transform(p).B = B;
   can_corrs.cc(p,:)            = r;
   can_corrs.lin_transform(p).U = U;
   can_corrs.lin_transform(p).V = V;
   can_corrs.stats(p)           = stats;
   if within_ccs
        can_corrs.cc_within_avg = cc_within_avg;
   end
end


%% ------------------------------------------------------------------------
% Summary stats

% DIMENSIONS FOR SUMMARY STATS
dims_stats      = 1:2;

diff_in_days    = days_btw_sessions(comb_tds);
diff_in_days    = diff_in_days(:,2) - diff_in_days(:,1);
 
mn_cc           = zeros(size(comb_tds,1),1);
sem_cc          = zeros(size(comb_tds,1),1);
mn_corr         = zeros(size(comb_tds,1),1);
sem_corr        = zeros(size(comb_tds,1),1);


for i = 1:size(comb_tds,1)
    cc          = abs(can_corrs.cc(i,dims_stats));
    cr          = abs(r_ctrl(i,dims_stats));
    mn_cc(i)    = mean(cc);
    sem_cc(i)   = std(cc)/sqrt(numel(dims));
    mn_corr(i)  = mean(cr);
    sem_corr(i) = std(cr)/sqrt(numel(dims));
end

% normalize corrs
if within_ccs
    mn_cc_norm      = zeros(size(comb_tds,1),1);
    sem_cc_norm     = zeros(size(comb_tds,1),1);
    mn_corr_norm    = zeros(size(comb_tds,1),1);
    sem_corr_norm   = zeros(size(comb_tds,1),1);


    for i = 1:size(comb_tds,1)
        cc_norm     = abs(can_corrs.cc(i,dims_stats))./can_corrs.cc_within_avg;
        cr_norm     = abs(r_ctrl(i,dims_stats))./can_corrs.cc_within_avg;
        mn_cc_norm(i) = mean(cc_norm);
        sem_cc_norm(i) = std(cc_norm)/sqrt(numel(dims));
        mn_corr_norm(i) = mean(cr_norm);
        sem_corr_norm(i) = std(cr_norm)/sqrt(numel(dims));
    end
end


%% ------------------------------------------------------------------------

% STANDARD CCA PLOT

% create legend
clear lg;
for i = 1:size(comb_tds,1)
    lg{i}           = [num2str(days_btw_sessions(comb_tds(i,1))) ' - ', ...
                        num2str(days_btw_sessions(comb_tds(i,2)))];
end

% plot!
cs = parula(size(comb_tds,1));
figure, hold on
for i = 1:size(comb_tds,1)
    plot(can_corrs.cc(i,:)','linewidth',2,'color',cs(i,:));
end
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Manifold projection')
ylabel('Canonical correlation')
legend(lg,'Location','NorthEast'), legend boxoff
xlim([1 max(dims)+3])
ylim([0 1])

% and the controls
for i = 1:size(comb_tds,1)
%    plot(abs(r_ctrl(i,:)'),'linewidth',2,'color',cs(i,:),'linestyle',':');
    plot(r_ctrl(i,:)','linewidth',2,'color',cs(i,:),'linestyle',':');
end


%% ------------------------------------------------------------------------

% BAR PLOT SUMMARIZING ACROSS DIMENSIONS

figure, hold on
for i = 1:size(comb_tds,1)
    errorbar((2*i-1),mn_cc(i),sem_cc(i),'k','linewidth',2)
    bar((2*i-1),mn_cc(i),'facecolor',cs(i,:));
    
    errorbar(2*i,mn_corr(i),sem_corr(i),'k','linewidth',2)
    bar(2*i,mn_corr(i),'facecolor',[1 1 1],'edgecolor',cs(i,:),'linewidth',4);
end
xlim([0 size(comb_tds,1)*2+1]),ylim([0 1]),title('Aligned (full) vs Not aligned (empty)')
set(gca,'TickDir','out','FontSize',14), box off, ylabel('Correlation')
set(gca,'XTick',1.5:2:size(comb_tds,1)*2,'XTickLabel',lg,'XTickLabelRotation',45)


% -------------------------------------------------------------------------
% Correlation as function of the time between sessions

figure, hold on
for i = 1:size(comb_tds,1)
    errorbar(diff_in_days(i),mn_cc(i),sem_cc(i),'.','markersize',36,'color',cs(i,:));
    errorbar(diff_in_days(i),mn_corr(i),sem_corr(i),'o','markersize',12,'color',cs(i,:));
end
set(gca,'TickDir','out','FontSize',14), box off, 
xlabel('Days between sessions'), ylabel('Correlation')
ylim([0 1]),xlim([0 max(diff_in_days)+1])
title('Aligned (full) vs Not aligned (empty)')


%% ------------------------------------------------------------------------

% SAME PLOTS WITH NORMALIZED CORRELATIONS

if within_ccs
    figure, hold on
    for i = 1:size(comb_tds,1)
        errorbar((2*i-1),mn_cc_norm(i),sem_cc_norm(i),'k','linewidth',2)
        bar((2*i-1),mn_cc_norm(i),'facecolor',cs(i,:));

        errorbar(2*i,mn_corr_norm(i),sem_corr_norm(i),'k','linewidth',2)
        bar(2*i,mn_corr_norm(i),'facecolor',[1 1 1],'edgecolor',cs(i,:),'linewidth',4);
    end
    xlim([0 size(comb_tds,1)*2+1]),ylim([0 1]),title('Aligned (full) vs Not aligned (empty)')
    set(gca,'TickDir','out','FontSize',14), box off, ylabel('Correlation')
    set(gca,'XTick',1.5:2:size(comb_tds,1)*2,'XTickLabel',lg,'XTickLabelRotation',45)


    % -------------------------------------------------------------------------
    % Correlation as function of the time between sessions

    figure, hold on
    for i = 1:size(comb_tds,1)
        errorbar(diff_in_days(i),mn_cc_norm(i),sem_cc_norm(i),'.','markersize',36,'color',cs(i,:));
        errorbar(diff_in_days(i),mn_corr_norm(i),sem_corr_norm(i),'o','markersize',12,'color',cs(i,:));
    end
    set(gca,'TickDir','out','FontSize',14), box off, 
    xlabel('Days between sessions'), ylabel('Correlation')
    ylim([0 1]),xlim([0 max(diff_in_days)+1])
    title('Aligned (full) vs Not aligned (empty)')
end