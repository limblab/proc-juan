%
% Figure 3 Aligned Latent activity
%


function SOT_Fig_3_aligned_latent_activity( td, align_results, meta, params, varargin )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
top_lv_plot             = 4;

if nargin > 1, assignParams(who, params); end % overwrite defaults

if nargin == 5
    within_day_align_results = varargin{1};
    mean_or_max_within_day  = 'max';
end

n_sessions = length(meta.sessions);



% -------------------------------------------------------------------------
% PLOT CC LEADING X LATENT VARIABLES VS NUMBER OF SESSIONS

latent_vars_plot = 1:top_lv_plot;
errorbars = 'sem'; % 'sd' 'sem'


dd = align_results.diff_days; % define for convenience


% Get some basic stats
mn_aligned_plot = arrayfun( @(x) mean(x.cc(latent_vars_plot)), align_results.aligned_info );
sd_aligned_plot = arrayfun( @(x) std(x.cc(latent_vars_plot)), align_results.aligned_info );
sem_aligned_plot = sd_aligned_plot/sqrt(length(latent_vars_plot));


mn_corr_plot = arrayfun( @(x) mean(abs(x.r(latent_vars_plot))), align_results.corr_info );
sd_corr_plot = arrayfun( @(x) std(abs(x.r(latent_vars_plot))), align_results.corr_info );
sem_corr_plot = sd_corr_plot/sqrt(length(latent_vars_plot));


% Linear fits to mean CCs / corrs
linfit_aligned = polyfit( dd, mn_aligned_plot, 1 );
linfit_unal = polyfit( dd, mn_corr_plot, 1 );
xfit = [0 max(dd)+1];
y_aligned = polyval(linfit_aligned,xfit);
y_unal = polyval(linfit_unal,xfit);


% retrieve data for within-day ceiling
if nargin == 5
    
    ceil_m = zeros(1,size(align_results.comb_sessions,1));
    ceil_eb = zeros(1,size(align_results.comb_sessions,1));
   
    cs = align_results.comb_sessions;
    
    warning('Not including errorbars for within-day latent activity alignment')
    
    % populate
    for c = 1:size(cs,1)
        
        s1 = cs(c,1);
        s2 = cs(c,2);
        
        switch mean_or_max_within_day
            case 'mean'
                aux_ceil = mean([within_day_align_results(s1).cc_m; within_day_align_results(s2).cc_m],1);
            case 'max'
                aux_ceil = max([within_day_align_results(s1).cc_m; within_day_align_results(s2).cc_m]);
        end
        
        ceil_m(c) = mean(aux_ceil(1:top_lv_plot));
        
        % ToDo
        switch errorbars
            case 'sd'
            case 'sem'
        end
    end
    
    % linear fit
    linfit_ceil = polyfit( dd, ceil_m, 1);
    y_ceil = polyval( linfit_ceil, xfit );
end


% distribution of stability over time
x_hist = 0:0.05:1.05;
n_hist = length(ceil_m);
y_hist_ceil = histcounts(ceil_m,x_hist)/n_hist*100;
y_hist_align = histcounts(mn_aligned_plot,x_hist)/n_hist*100;
y_hist_corr = histcounts(mn_corr_plot,x_hist)/n_hist*100;

mn_hist_ceil = mean(ceil_m);
mn_hist_align = mean(mn_aligned_plot);
mn_hist_corr = mean(mn_corr_plot);

sd_hist_ceil = std(ceil_m);
sd_hist_align = std(mn_aligned_plot);
sd_hist_corr = std(mn_corr_plot);

% sem_hist_ceil = std(ceil_m)/sqrt(numel(ceil_m));
% sem_hist_align = std(mn_aligned_plot)/sqrt(numel(mn_aligned_plot));
% sem_hist_corr = std(mn_corr_plot)/sqrt(numel(mn_corr_plot));

ymax_hist = max([y_hist_ceil, y_hist_align, y_hist_corr]');

y2_hist_fig = ceil(ymax_hist/10)*10;
ystats_hist = (y2_hist_fig - ymax_hist)/2 + ymax_hist; 

% NON-parametric comparison unaligned and aligned
p_align_within = ranksum(ceil_m,mn_aligned_plot);
p_align_unaligned = ranksum(mn_corr_plot,mn_aligned_plot);


% The same for the normalized aligned and unaligned data
norm_align = mn_aligned_plot./ceil_m;
norm_corr = mn_corr_plot./ceil_m;
y_hist_norm_align = histcounts(norm_align,x_hist)/n_hist*100;
y_hist_norm_corr = histcounts(norm_corr,x_hist)/n_hist*100;

mn_norm_align = mean(norm_align);
mn_norm_corr = mean(norm_corr);
sd_norm_align = std(norm_align);
sd_norm_corr = std(norm_corr);

y_max_nh = max([y_hist_norm_align,y_hist_norm_corr]);
y2_normhist = ceil((y_max_nh)/10)*10+5;
ystats_normhist = (y2_normhist - y_max_nh)/2 + y_max_nh;

% ---------------------------------------------------------------------------------
% Plots vs time

col_unal = [.6 .6 .6];
col_ceil = [1 0.6 0];


f1 = figure; hold on
plot(xfit,y_aligned,'k','linewidth',2)
plot(xfit,y_unal,'color',col_unal,'linewidth',2)
if exist('within_day_align_results','var')
    plot(xfit,y_ceil,'color',col_ceil,'linewidth',2)
    plot(dd,ceil_m,'linestyle','none','marker','.','markersize',32,'color',col_ceil)
end
switch errorbars
    case 'sd'
        errorbar(dd,mn_aligned_plot,sd_aligned_plot,'k','linewidth',1.5,'linestyle','none','marker','.','markersize',32)
        errorbar(dd,mn_corr_plot,sd_corr_plot,'color',col_unal,'linewidth',1.5,'linestyle','none','marker','.','markersize',32)
    case 'sem'
        errorbar(dd,mn_aligned_plot,sem_aligned_plot,'k','linewidth',1.5,'linestyle','none','marker','.','markersize',32)
        errorbar(dd,mn_corr_plot,sem_corr_plot,'color',col_unal,'linewidth',1.5,'linestyle','none','marker','.','markersize',32)
end
set(gca,'TickDir','out','FontSize',14), box off
set(gcf, 'color', [1 1 1])
if exist('within_day_align_results','var')
    legend('Aligned','Unaligned','Within-day'), 
else
    legend('Aligned','Unaligned'), 
end
legend boxoff
ylim([0 1]), xlabel('Days from first session'), ylabel(['CC latent activity (top ' num2str(length(latent_vars_plot)) ' modes)'])







% % -------------------------------------------------------------------------
% % PLOT HISTOGRAM ALIGNED NOT ALIGNED CCs
% 
% bin_hist = 0.05;
% x_hist = 0:bin_hist:(1+bin_hist);
% 
% data_aligned = cell2mat(arrayfun( @(x) x.cc(latent_vars_plot), align_results.aligned_info, 'UniformOutput', false));
% data_unaligned = cell2mat(arrayfun( @(x) abs(x.r(latent_vars_plot)), align_results.corr_info, 'UniformOutput', false));
% 
% hist_aligned = histcounts(data_aligned,x_hist)/length(data_aligned)*100;
% hist_unaligned = histcounts(data_unaligned,x_hist)/length(data_unaligned)*100;
% 
% figure, hold on
% hha = bar(x_hist(1:end-1),hist_aligned,'histc');
% set(hha,'facecolor','k')
% hhu = bar(x_hist(1:end-1),hist_unaligned,'histc');
% set(hhu,'facecolor',col_unal);
% alpha(hhu,0.5)
% set(gca,'TickDir','out','FontSize',14), box off
% set(gcf, 'color', [1 1 1])
% legend('Aligned','Unaligned','Location','NorthWest'), legend boxoff
% xlim([0 1]), xlabel('CC (aligned) or corr (unaligned) across days')
% ylabel('Mode comparisons (%)')



%% % ----------------------------------------------------------------------
% PLOT THE MEAN ALIGNED MEAN TRAJECTORIES PER TARGET

s1 = 1;
s2 = n_sessions;

diff_days = datenum(meta.sessions{s2})-datenum(meta.sessions{s1});

modesplot = 1:3;


% Get the comparison that we want from align_results
idx_cmp = find( align_results.comb_sessions(:,1)==s1 & align_results.comb_sessions(:,2)==s2 );

% Hack: create aux TD signals with each of these sessions
[~, td1] = getTDidx( td, {'date',meta.sessions(s1)} );
[~, td2] = getTDidx( td, {'date',meta.sessions(s2)} );

% define name var to plot
latent_signals_plot = align_results.signals;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate concatenated aligned latent activity into trials and
% add them to TD
bins_p_trial = size(td1(1).pos,1);

for t = 1:length(td1)
    be = bins_p_trial*(t-1)+1;
    en = bins_p_trial*t;
    td1(t).(latent_signals_plot) = align_results.aligned_info(idx_cmp).U(be:en,:);
    td2(t).(latent_signals_plot) = align_results.aligned_info(idx_cmp).V(be:en,:);
end



% For the aligned data --stored in struct td2
td_avg_ali1  = trialAverage(td1,'target_direction');
td_avg_ali2  = trialAverage(td2,'target_direction');

 
cols        = parula(length(td_avg_ali1));
 
f3 = figure; hold on
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
set(f3,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title('Day 1 Aligned')
ax1 = gca;


f4 = figure; hold on
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
set(f4,'color', [1 1 1]);
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

figure(f3), xlim(xl); ylim(yl); zlim(zl);
figure(f4), xlim(xl); ylim(yl); zlim(zl);



% HISTOGRAMS

f5 = figure; hold on
b1 = bar(x_hist(1:end-1),y_hist_corr,'histc');
set(b1,'facecolor',col_unal)
alpha(b1,0.5)
b2 = bar(x_hist(1:end-1),y_hist_ceil,'histc');
set(b2,'facecolor',col_ceil)
alpha(b2,0.5)
b3 = bar(x_hist(1:end-1),y_hist_align,'histc');
set(b3,'facecolor','k')
alpha(b3,0.5)
set(gca,'TickDir','out','FontSize',14), box off
set(gcf, 'color', [1 1 1])
xlabel(['CC latent activity (top ' num2str(length(latent_vars_plot)) ' modes)'])
ylabel('Session comparisons (%)')
xlim([0 1]); ylim([0 y2_hist_fig])
plot(mn_hist_ceil,ystats_hist,'.','color',col_ceil,'markersize',32);
plot([mn_hist_ceil-sd_hist_ceil, mn_hist_ceil+sd_hist_ceil],[ystats_hist ystats_hist],'color',col_ceil,'linewidth',2);
plot(mn_hist_corr,ystats_hist,'.','color',col_unal,'markersize',32);
plot([mn_hist_corr-sd_hist_corr, mn_hist_corr+sd_hist_corr],[ystats_hist ystats_hist],'color',col_unal,'linewidth',2);
plot(mn_hist_align,ystats_hist,'.','color','k','markersize',32);
plot([mn_hist_align-sd_hist_align, mn_hist_align+sd_hist_align],[ystats_hist ystats_hist],'color','k','linewidth',2);
legend('Unaligned','Within-day','Aligned','Location','West'),legend boxoff

text(.1,ystats_hist-5,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,ystats_hist-10,['P=' num2str(p_align_unaligned,2)],'Fontsize',14)


% Normalized hist
f6 = figure; hold on
bn1 = bar(x_hist(1:end-1),y_hist_norm_corr,'histc');
set(bn1,'facecolor',col_unal)
alpha(bn1,0.5)
bn2 = bar(x_hist(1:end-1),y_hist_norm_align,'histc');
set(bn2,'facecolor','k')
alpha(bn2,0.5)
set(gca,'TickDir','out','FontSize',14), box off
set(gcf, 'color', [1 1 1])
xlabel(['Normalized CC latent activity (top ' num2str(length(latent_vars_plot)) ' modes)'])
ylabel('Session comparisons (%)')
xlim([0 1]); ylim([0 y2_normhist])
plot(mn_norm_corr,ystats_normhist,'.','color',col_unal,'markersize',32);
plot([mn_norm_corr-sd_norm_corr, mn_norm_corr+sd_norm_corr],[ystats_normhist ystats_normhist],'color',col_unal,'linewidth',2);
plot(mn_norm_align,ystats_normhist,'.','color','k','markersize',32);
plot([mn_norm_align-sd_norm_align, mn_norm_align+sd_norm_align],[ystats_normhist ystats_normhist],'color','k','linewidth',2);
legend('Within-day','Aligned','Location','West'),legend boxoff

text(.1,ystats_hist-5,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,ystats_hist-8,['P=' num2str(p_align_unaligned,2)],'Fontsize',14)





% SAVE FIG?
if params.save_fig
    
    % Stability over time
    fn1 = [td(1).monkey '_' params.signals(1:end-4) '_Latent_activity_over_time_' num2str(length(params.mani_dims)) 'D'];

    savefig(f1,fullfile(params.save_dir,params.signals(1:end-4),[fn1 '.fig']));
    saveas(f1,fullfile(params.save_dir,params.signals(1:end-4),[fn1 '.png']));
    saveas(f1,fullfile(params.save_dir,params.signals(1:end-4),[fn1 '.pdf']));
    
    
    % Normalized Histograms
    fn2 = [td(1).monkey '_' params.signals(1:end-4) '_Latent_activity_normalized_distribution_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f6,fullfile(params.save_dir,params.signals(1:end-4),[fn2 '.fig']));
    saveas(f6,fullfile(params.save_dir,params.signals(1:end-4),[fn2 '.png']));
    saveas(f6,fullfile(params.save_dir,params.signals(1:end-4),[fn2 '.pdf']));
    
    
    % Aligned latent trajectories
    fn3 = [td(1).monkey '_' params.signals(1:end-4) '_Aligned_latent_trajectories_Day_1_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f4,fullfile(params.save_dir,params.signals(1:end-4),[fn3 '.fig']));
    saveas(f4,fullfile(params.save_dir,params.signals(1:end-4),[fn3 '.png']));
    saveas(f4,fullfile(params.save_dir,params.signals(1:end-4),[fn3 '.pdf']));
    
    fn4 = [td(1).monkey '_' params.signals(1:end-4) '_Aligned_latent_trajectories_Day_' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f5,fullfile(params.save_dir,params.signals(1:end-4),[fn4 '.fig']));
    saveas(f5,fullfile(params.save_dir,params.signals(1:end-4),[fn4 '.png']));
    saveas(f5,fullfile(params.save_dir,params.signals(1:end-4),[fn4 '.pdf']));
    
end