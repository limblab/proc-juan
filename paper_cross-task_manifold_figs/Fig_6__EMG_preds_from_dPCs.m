%
% Figure 6: EMG predictions from dPCs. 
%
%


% Need to call comp_decoders_dPCA.m first 
if ~exist('dPCA_fit','var')
    comp_decoders_dPCA;
end

% close all the previously generated figures
% close all;

ds_ex       = 4;
t_win_ex    = [5 35];
ex_emg      = [4 6]; % Among the chosen ones
% margs: {'task','target','time','task/tgt int'};
marg_ex     = 2;


% -------------------------------------------------------------------------
% Some preliminary stuff for the predictions


% create time vectors
bin_size    = round(mean(diff(datasets{1}.binned_data{1}.timeframe))*100)/100;

t           = 0 : bin_size : bin_size * (size(dPCA_fit{ds_ex,marg_ex}.all{1}.y,1)-1);

v_lines_idx = dPCA_fit{ds_ex,marg_ex}.trial_beg-5;
v_lines_idx(1) = [];

% get indexes for plot window
idx_1       = find( t >= t_win_ex(1),1 );
idx_2       = find( t >= t_win_ex(2),1 );


% get raw EMGs and predictions
y_1         = dPCA_fit{ds_ex,marg_ex}.all{1}.y(:,emg_chs{ds_ex}(ex_emg(1)));
y_2         = dPCA_fit{ds_ex,marg_ex}.all{1}.y(:,emg_chs{ds_ex}(ex_emg(2)));

pred_1      = dPCA_fit{ds_ex,marg_ex}.all{1}.y_pred(:,emg_chs{ds_ex}(ex_emg(1)));
pred_2      = dPCA_fit{ds_ex,marg_ex}.all{1}.y_pred(:,emg_chs{ds_ex}(ex_emg(2)));


% Compute R2's in analysis window for the plot that shows raw predictions
R2_win_1    = CalculateR2( y_1(idx_1:idx_2), pred_1(idx_1:idx_2) );
R2_win_2    = CalculateR2( y_2(idx_1:idx_2), pred_2(idx_1:idx_2) );

R2_tot_1    = dPCA_fit{ds_ex,marg_ex}.all{1}.R2(emg_chs{ds_ex}(ex_emg(1)));
R2_tot_2    = dPCA_fit{ds_ex,marg_ex}.all{1}.R2(emg_chs{ds_ex}(ex_emg(2)));



% -------------------------------------------------------------------------
% Create data matrices for distributions for boxplots

norm_R2_marg_all_mtrx = cell2mat(norm_R2_marg_all')';



% -------------------------------------------------------------------------
% The figure

% top panel (1) example predictions
figure('units','normalized','outerposition',[0.25 0.3 0.5 0.4])
hold on
plot(t,y_1,'color',[.5 .5 .5],'linewidth',1)
plot(t,pred_1,'k','linewidth',1)
legend({['R^2_{win}=' num2str(R2_win_1,2)],['R^2_{tot}=' num2str(R2_tot_1,2)]})
xlim(t_win_ex),ylim([0 1.2])
for v = 1:length(v_lines_idx)
    plot([t(v_lines_idx(v)),t(v_lines_idx(v))],[0 1.5],'color',[1 1 1],'linewidth',3)%[.5 .5 .5],'linewidth',1)
end 
set(gca,'TickDir','out','FontSize',12), box off
ylabel(datasets{D(ds_ex)}.binned_data{1}.emgguide(emg_chs{ds_ex}(ex_emg(1))))   

% top panel (2) another example pred
figure('units','normalized','outerposition',[0.25 0.3 0.5 0.4])
hold on
plot(t,y_2,'color',[.5 .5 .5],'linewidth',1)
plot(t,pred_2,'k','linewidth',1)
legend({['R^2_{win}=' num2str(R2_win_2,2)],['R^2_{tot}=' num2str(R2_tot_2,2)]})
xlim(t_win_ex),ylim([0 1.2])
for v = 1:length(v_lines_idx)
    plot([t(v_lines_idx(v)),t(v_lines_idx(v))],[0 1.5],'color',[1 1 1],'linewidth',3),%[.5 .5 .5],'linewidth',1)
end   
set(gca,'TickDir','out','FontSize',12), box off
ylabel(datasets{D(ds_ex)}.binned_data{1}.emgguide(emg_chs{ds_ex}(ex_emg(2))))
xlabel('Time (s)')


% ORIGINAL VERSION

figure('units','normalized','outerposition',[0.25 0.25 0.5 0.45])
% bottom left panel: Norm R^2 preds
subplot(1,2,1), hold on
for s = 1:4
    if use_sem_yn, eb = sem_norm_R2_marg_all(s); 
        else eb = sd_norm_R2_marg_all(s); end
    % and plot
    errorbar(s,mn_norm_R2_marg_all(s), eb,...
        'marker','none','color',marg_cols(s,:),...
        'linewidth',1,'linestyle','none')
    bar(s,mn_norm_R2_marg_all(s),'FaceColor',marg_cols(s,:))
end
ylim([0 1]); xlim([0 5]); box off
ylabel('Norm. R^2 EMG')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:s,'XTickLabel',{'Task','Target','Dynamics','Task-target'},'XTickLabelRotation',45)

% bottom right panel: Norm R^2 preds vs. Neural var explained
subplot(1,2,2), hold on
for k = 1:numel(monk)
    for m = 1:4
        if k == 1 % Jaco
            plot(neural_marg_mn{k}(m),mn_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor','none')
        else
            plot(neural_marg_mn{k}(m),mn_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor',marg_cols(m,:))
        end
        plot([neural_marg_mn{k}(m)-neural_marg_sd{k}(m),...
                neural_marg_mn{k}(m)+neural_marg_sd{k}(m)],...
                [mn_norm_R2_marg_all_musc_per_monk(k,m), mn_norm_R2_marg_all_musc_per_monk(k,m)],...
                'color',marg_cols(m,:))
        plot([neural_marg_mn{k}(m), neural_marg_mn{k}(m)],...
            [mn_norm_R2_marg_all_musc_per_monk(k,m)-sd_norm_R2_marg_all_musc_per_monk(k,m),...
            mn_norm_R2_marg_all_musc_per_monk(k,m)+sd_norm_R2_marg_all_musc_per_monk(k,m)],...
            'color',marg_cols(m,:))
    end
end
ylim([0 1]);xlim([0 .45])
ylabel('R^2 EMG')
xlabel('Neural varance expl. (%)')
set(gca,'TickDir','out','FontSize',12)


% UPDATED VERSION WITH BOX PLOT
figure('units','normalized','outerposition',[0.25 0.25 0.5 0.45])
% bottom left panel: Norm R^2 preds
subplot(1,2,1), hold on
for s = 1:4
    % boxplot(norm_R2_marg_all_mtrx); %,'color',marg_cols(s,:));
    bplot(norm_R2_marg_all_mtrx(:,s),s,'nomean')%,'color',marg_cols(s,:));
end
ylim([0 1]); xlim([0 5]); box off
ylabel('Norm. R^2 EMG')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:s,'XTickLabel',{'Task','Target','Dynamics','Task-target'},'XTickLabelRotation',45)

% bottom right panel: Norm R^2 preds vs. Neural var explained
subplot(1,2,2), hold on
for k = 1:numel(monk)
    for m = 1:4
        if k == 1 % Jaco
            plot(neural_marg_median{k}(m),median_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor','none')
        else
            plot(neural_marg_median{k}(m),median_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor',marg_cols(m,:))
        end
        plot([neural_marg_median{k}(m)-neural_marg_sd{k}(m),...
                neural_marg_median{k}(m)+neural_marg_sd{k}(m)],...
                [median_norm_R2_marg_all_musc_per_monk(k,m), median_norm_R2_marg_all_musc_per_monk(k,m)],...
                'color',marg_cols(m,:))
        plot([neural_marg_median{k}(m), neural_marg_median{k}(m)],...
            [median_norm_R2_marg_all_musc_per_monk(k,m)-sd_norm_R2_marg_all_musc_per_monk(k,m),...
            median_norm_R2_marg_all_musc_per_monk(k,m)+sd_norm_R2_marg_all_musc_per_monk(k,m)],...
            'color',marg_cols(m,:))
    end
end
ylim([0 1]);xlim([0 .45])
ylabel('R^2 EMG')
xlabel('Neural varance expl. (%)')
set(gca,'TickDir','out','FontSize',12)
