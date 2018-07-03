%
% Suppl Fig 9d: Predictions of EMG dPCs based on different sets of neural dPCs
% -- Need to run predict_EMG_dPCs_from_neural_dPCs
%
% TODO: Need to double check how far the whiskers extend to
%


% paper colors
marg_cols       = [.4 .6 .5; .7 .1 .1; .5 .65 .9; 1 .6 .3];

% color individual datapoints
pt_cols         = [.75 .75 .75];

% marg labels
marg_labels     = {'task','target','time','task-target'};


figure, hold on
%boxplot(norm_R2_marg);
% for m = 1:size(norm_R2_marg,2)
%     bplot(norm_R2_marg(~isnan(norm_R2_marg(:,m)),m),m,'nomean')
% end
boxplot(norm_R2_marg,'Whisker',.5)
ylim([0 1]); xlim([0 5]); box off
ylabel('Norm. R^2 EMG')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:m,'XTickLabel',{'Task','Target','Dynamics','Task-target'},'XTickLabelRotation',45)

