%
% Plot "stability" of the behavior
%

function SOT_Fig_stability_behavior( corr_kin, params )


% Linear fits
x_lf    = [0 max(corr_kin.diff_days)+1];

lfx     = polyfit(corr_kin.diff_days,corr_kin.r(:,1)',1);
lfy     = polyfit(corr_kin.diff_days,corr_kin.r(:,2)',1);

y_lfx   = polyval(lfx,x_lf);
y_lfy   = polyval(lfy,x_lf);

cols    = parula(3);


% distribution of mean X-Y corr
m_r     = mean(corr_kin.r,2);
x_hist  = 0:0.02:1.02;
y_hist  = histcounts(m_r,x_hist)/numel(m_r)*100;
mn_hist = mean(m_r);
sem_hist = std(m_r)/sqrt(length(m_r));
sd_hist = std(m_r);

yl_hist = [0 floor((max(y_hist)+10))];
y_stats = floor((max(y_hist)+5));


% PLOT CORR VS TIME
f1 = figure; hold on;
plot(x_lf,y_lfx,'linewidth',2,'color',cols(1,:))
plot(x_lf,y_lfy,'linewidth',2,'color',cols(2,:))
plot(corr_kin.diff_days,corr_kin.r(:,1),'.','markersize',32,'color',cols(1,:))
plot(corr_kin.diff_days,corr_kin.r(:,2),'.','markersize',32,'color',cols(2,:))
ylim([0 1])
set(gca,'TickDir','out','FontSize',14), box off
legend(['X ' corr_kin.var],['Y ' corr_kin.var],'Location','SouthEast'),legend boxoff 
xlabel('Days between sessions'),
if ~params.stab_behav.trial_avg
    ylabel(['Correlation hand ' corr_kin.var])
else
    ylabel(['Correlation hand ' corr_kin.var ' trial-averaged'])
end
set(gcf,'color','w')


% PLOT DISTRIBUTION MEAN CORRS
f2 = figure; hold on;
hb = bar(x_hist(1:end-1),y_hist,'histc');
plot(mn_hist,y_stats,'.','color',[.7 .7 .7],'markersize',32);
plot([mn_hist+sd_hist mn_hist-sd_hist],[y_stats y_stats],'color',[.7 .7 .7],'linewidth',2);
% plot(sem_hist,y_stats,'color',[.7 .7 .7],'linewidth',2);
set(hb,'FaceColor',[.7 .7 .7]);
set(gca,'TickDir','out','FontSize',14), box off
ylim(yl_hist); xlim([0 1]);
if ~params.stab_behav.trial_avg
    xlabel(['Mean correlation X-Y hand ' corr_kin.var])
else
    xlabel(['Mean correlation X-Y hand ' corr_kin.var ' trial-averaged'])
end
ylabel('Session comparisons (%)')
set(gcf,'color','w')


% SAVE FIG ?
if params.stab_behav.save_fig
    
    ff = '/Users/juangallego/Dropbox/Juan and Matt - Stability latent activity/Results/Behavior';
    fn1 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Behavior_over_time'];

    savefig(f1,[ff filesep fn1]);
    saveas(f1,[ff filesep fn1 '.png']);
    saveas(f1,[ff filesep fn1 '.pdf']);
    
    fn2 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Behavior_distribution'];
    
    savefig(f2,[ff filesep fn2]);
    saveas(f2,[ff filesep fn2 '.png']);
    saveas(f2,[ff filesep fn2 '.pdf']);
end
