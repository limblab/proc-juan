%
% Plot "stability" of the behavior
%


x_lf = [0 max(corr_kin.diff_days)+1];

lfx = polyfit(corr_kin.diff_days,corr_kin.r(:,1)',1);
lfy = polyfit(corr_kin.diff_days,corr_kin.r(:,2)',1);

y_lfx = polyval(lfx,x_lf);
y_lfy = polyval(lfy,x_lf);

cols = parula(3);

figure,hold on
plot(x_lf,y_lfx,'linewidth',2,'color',cols(1,:))
plot(x_lf,y_lfy,'linewidth',2,'color',cols(2,:))
plot(corr_kin.diff_days,corr_kin.r(:,1),'.','markersize',32,'color',cols(1,:))
plot(corr_kin.diff_days,corr_kin.r(:,2),'.','markersize',32,'color',cols(2,:))
ylim([0 1])
set(gca,'TickDir','out','FontSize',14), box off
legend(['X ' corr_kin.var],['Y ' corr_kin.var],'Location','SouthEast'),legend boxoff 
xlabel('Days between sessions'),ylabel(['Corr hand ' corr_kin.var])
set(gcf,'color','w')
