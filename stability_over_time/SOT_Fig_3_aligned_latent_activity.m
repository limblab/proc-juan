%
% Figure 3 Aligned Latent activity
%


% -------------------------------------------------------------------------
% PLOT CC LEADING X LATENT VARIABLES VS NUMBER OF SESSIONS

latent_vars_plot = 1:4;
errorbars = 'sem'; % 'sd' 'sem'


dd = aligned_latent_results.diff_days; % define for convenience


% Get some basic stats
mn_aligned_plot = arrayfun( @(x) mean(x.cc(latent_vars_plot)), aligned_latent_results.aligned_info );
sd_aligned_plot = arrayfun( @(x) std(x.cc(latent_vars_plot)), aligned_latent_results.aligned_info );
sem_aligned_plot = sd_aligned_plot/sqrt(length(latent_vars_plot));


mn_corr_plot = arrayfun( @(x) mean(abs(x.r(latent_vars_plot))), aligned_latent_results.corr_info );
sd_corr_plot = arrayfun( @(x) std(abs(x.r(latent_vars_plot))), aligned_latent_results.corr_info );
sem_corr_plot = sd_corr_plot/sqrt(length(latent_vars_plot));


% Linear fits to mean CCs / corrs
linfit_aligned = polyfit( dd, mn_aligned_plot, 1 );
linfit_unal = polyfit( dd, mn_corr_plot, 1 );
xfit = [0 max(dd)+1];
y_aligned = polyval(linfit_aligned,xfit);
y_unal = polyval(linfit_unal,xfit);


% Plot
figure, hold on
plot(xfit,y_aligned,'k','linewidth',2)
plot(xfit,y_unal,'color',[.65 .65 .65],'linewidth',2)
switch errorbars
    case 'sd'
        errorbar(dd,mn_aligned_plot,sd_aligned_plot,'k','linewidth',1.5,'linestyle','none','marker','.','markersize',32)
        errorbar(dd,mn_corr_plot,sd_corr_plot,'color',[.65 .65 .65],'linewidth',1.5,'linestyle','none','marker','.','markersize',32)
    case 'sem'
        errorbar(dd,mn_aligned_plot,sem_aligned_plot,'k','linewidth',1.5,'linestyle','none','marker','.','markersize',32)
        errorbar(dd,mn_corr_plot,sem_corr_plot,'color',[.65 .65 .65],'linewidth',1.5,'linestyle','none','marker','.','markersize',32)
end
set(gca,'TickDir','out','FontSize',14), box off
set(gcf, 'color', [1 1 1])
legend('Aligned','Unaligned'), legend boxoff
ylim([0 1]), xlabel('Days from first session'), ylabel(['CC latent activity (top ' num2str(length(latent_vars_plot)) ' modes)'])




% % -------------------------------------------------------------------------
% % PLOT HISTOGRAM ALIGNED NOT ALIGNED CCs
% 
% bin_hist = 0.05;
% x_hist = 0:bin_hist:(1+bin_hist);
% 
% data_aligned = cell2mat(arrayfun( @(x) x.cc(latent_vars_plot), aligned_latent_results.aligned_info, 'UniformOutput', false));
% data_unaligned = cell2mat(arrayfun( @(x) abs(x.r(latent_vars_plot)), aligned_latent_results.corr_info, 'UniformOutput', false));
% 
% hist_aligned = histcounts(data_aligned,x_hist)/length(data_aligned)*100;
% hist_unaligned = histcounts(data_unaligned,x_hist)/length(data_unaligned)*100;
% 
% figure, hold on
% hha = bar(x_hist(1:end-1),hist_aligned,'histc');
% set(hha,'facecolor','k')
% hhu = bar(x_hist(1:end-1),hist_unaligned,'histc');
% set(hhu,'facecolor',[.65 .65 .65]);
% alpha(hhu,0.5)
% set(gca,'TickDir','out','FontSize',14), box off
% set(gcf, 'color', [1 1 1])
% legend('Aligned','Unaligned','Location','NorthWest'), legend boxoff
% xlim([0 1]), xlabel('CC (aligned) or corr (unaligned) across days')
% ylabel('Mode comparisons (%)')



clear mn* sd* sem* errorbars y_* xfit linfit* latent_vars_plot ans



%% % ----------------------------------------------------------------------
% PLOT THE MEAN ALIGNED MEAN TRAJECTORIES PER TARGET

s1 = 1;
s2 = n_sessions;

diff_days = datenum(meta.sessions{s2})-datenum(meta.sessions{s1});

modesplot = 1:3;


% Get the comparison that we want from aligned_latent_results
idx_cmp = find( aligned_latent_results.comb_sessions(:,1)==s1 & aligned_latent_results.comb_sessions(:,2)==s2 );

% Hack: create aux TD signals with each of these sessions
[~, td1] = getTDidx( master_td, {'date',meta.sessions(s1)} );
[~, td2] = getTDidx( master_td, {'date',meta.sessions(s2)} );

% define name var to plot
latent_signals_plot = ['align_' pars.spiking_inputs{1}(1:end-7) '_pca'];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate concatenated aligned latent activity into trials and
% add them to TD
bins_p_trial = size(td1(1).pos,1);

for t = 1:length(td1)
    be = bins_p_trial*(t-1)+1;
    en = bins_p_trial*t;
    td1(t).(latent_signals_plot) = aligned_latent_results.aligned_info(idx_cmp).U(be:en,:);
    td2(t).(latent_signals_plot) = aligned_latent_results.aligned_info(idx_cmp).V(be:en,:);
end



% For the aligned data --stored in struct td2
td_avg_ali1  = trialAverage(td1,'target_direction');
td_avg_ali2  = trialAverage(td2,'target_direction');

 
cols        = parula(length(td_avg_ali1));
 
f1 = figure; hold on
for t = 1:length(cols)
    plot3(td_avg_ali1(t).(latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali1(t).(latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali1(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali1(t).(latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali1(t).(latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali1(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(-169,13)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f1,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title('Day 1 Aligned')
ax1 = gca;


f2 = figure; hold on
for t = 1:length(cols)
    plot3(td_avg_ali2(t).(latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali2(t).(latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali2(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali2(t).(latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali2(t).(latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali2(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(-169,13)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f2,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title(['Day ' num2str(diff_days) ' Aligned'])
ax2 = gca;

% Make axes lims equal
xl(1) = min(ax1.XLim(1),ax2.XLim(1));
yl(1) = min(ax1.YLim(1),ax2.YLim(1));
zl(1) = min(ax1.ZLim(1),ax2.ZLim(1));
    
xl(2) = min(ax1.XLim(2),ax2.XLim(2));
yl(2) = min(ax1.YLim(2),ax2.YLim(2));
zl(2) = min(ax1.ZLim(2),ax2.ZLim(2));

figure(f1), xlim(xl); ylim(yl); zlim(zl);
figure(f2), xlim(xl); ylim(yl); zlim(zl);


clear dd idx_cmp td1 td2 td_avg* cols modesplot diff_days xl yl zl f1 f2 ax* be en bins_p_trial s1 s2 t this_s ans;