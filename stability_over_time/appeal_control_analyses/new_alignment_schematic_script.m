function new_alignment_schematic_script( td, save_dir, align_results, meta, params, varargin )
close all;

top_lv_plot             = 4;
modesplot = 1:3;

assignParams(who, params);

if nargin == 6
    within_day_align_results = varargin{1};
    mean_or_max_within_day  = 'max';
end

latent_vars_plot = 1:top_lv_plot;



%% prepare for plotting
s1 = 5;
s2 = 14;%n_sessions;

diff_days = datenum(meta.sessions{s2})-datenum(meta.sessions{s1});


% Get the comparison that we want from align_results
idx_cmp = find( align_results.comb_sessions(:,1)==s1 & align_results.comb_sessions(:,2)==s2 );

cc = align_results.cc(idx_cmp,:);
r = align_results.r(idx_cmp,:);
A = align_results.aligned_info(idx_cmp).A;
B = align_results.aligned_info(idx_cmp).B;





% Hack: create aux TD signals with each of these sessions
[~, td1] = getTDidx( td, {'date',meta.sessions(s1)} );
[~, td2] = getTDidx( td, {'date',meta.sessions(s2)} );

% define name var to plot
latent_signals_plot = align_results.signals;
align_latent_signals_plot = [latent_signals_plot '_align'];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate concatenated aligned latent activity into trials and
% add them to TD
bins_p_trial = size(td1(1).pos,1);

for t = 1:length(td1)
    be = bins_p_trial*(t-1)+1;
    en = bins_p_trial*t;
    td1(t).(align_latent_signals_plot) = align_results.aligned_info(idx_cmp).U(be:en,:);
    td2(t).(align_latent_signals_plot) = align_results.aligned_info(idx_cmp).V(be:en,:);
end



% For the aligned data --stored in struct td2
td_avg_ali1  = trialAverage(td1,'target_direction');
td_avg_ali2  = trialAverage(td2,'target_direction');


cols        = parula(length(td_avg_ali1)+1);



%% plot aligned and unaligned trajectories
% v1 = 169;
% v2 = 13;
v1 = -109;
v2 = 45;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALIGNED TRAJECTORIES
f1 = figure; hold on
for t = 1:length(td_avg_ali1)
    plot3(td_avg_ali1(t).(align_latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali1(t).(align_latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali1(t).(align_latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali1(t).(align_latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali1(t).(align_latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali1(t).(align_latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(v1,v2)
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
for t = 1:length(td_avg_ali1)
    plot3(td_avg_ali2(t).(align_latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali2(t).(align_latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali2(t).(align_latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali2(t).(align_latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali2(t).(align_latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali2(t).(align_latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(v1,v2)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f2,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title(['Day ' num2str(diff_days) ' Aligned'])
ax2 = gca;

% UNALIGNED TRAJECTORIES
f3 = figure; hold on
for t = 1:length(td_avg_ali1)
    plot3(td_avg_ali1(t).(latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali1(t).(latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali1(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali1(t).(latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali1(t).(latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali1(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(v1,v2)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f3,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title('Day 1 Unaligned')
ax3 = gca;



f4 = figure; hold on
for t = 1:length(td_avg_ali1)
    plot3(td_avg_ali2(t).(latent_signals_plot)(:,modesplot(1)),...
        td_avg_ali2(t).(latent_signals_plot)(:,modesplot(2)),...
        td_avg_ali2(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot3(td_avg_ali2(t).(latent_signals_plot)(1,modesplot(1)),...
        td_avg_ali2(t).(latent_signals_plot)(1,modesplot(2)),...
        td_avg_ali2(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    view(v1,v2)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f4,'color', [1 1 1]);
xlabel(['Neural mode ' num2str(modesplot(1))]);
ylabel(['Neural mode ' num2str(modesplot(2))]);
zlabel(['Neural mode ' num2str(modesplot(3))]);
grid on
title(['Day ' num2str(diff_days) ' Unaligned'])
ax4 = gca;

% Make axes lims equal
xl(1) = min([ax1.XLim(1),ax2.XLim(1),ax3.XLim(1),ax4.XLim(1)]);
yl(1) = min([ax1.YLim(1),ax2.YLim(1),ax3.YLim(1),ax4.YLim(1)]);
zl(1) = min([ax1.ZLim(1),ax2.ZLim(1),ax3.ZLim(1),ax4.ZLim(1)]);

xl(2) = max([ax1.XLim(2),ax2.XLim(2),ax3.XLim(2),ax4.XLim(2)]);
yl(2) = max([ax1.YLim(2),ax2.YLim(2),ax3.YLim(2),ax4.YLim(2)]);
zl(2) = max([ax1.ZLim(2),ax2.ZLim(2),ax3.ZLim(2),ax4.ZLim(2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% now plot the time varying projections

% UNALIGNED TRAJECTORIES
f5 = figure('Position',[100 100 800 500]); hold on
td_use = td_avg_ali1;
for t = 1:length(td_use)
    subplot(3,4,1); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(1)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(1)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,5); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(2)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(2)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,9); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
end


% now do day 2
td_use = td_avg_ali2;
for t = 1:length(td_use)
    subplot(3,4,2); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(1)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(1)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,6); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(2)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(2)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,10); hold all;
    plot(td_use(t).(latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
end







% ALIGNED TRAJECTORIES
td_use = td_avg_ali1;
for t = 1:length(td_use)
    subplot(3,4,3); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(1)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(1)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,7); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(2)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(2)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,11); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
end

% now do day 2
td_use = td_avg_ali2;
for t = 1:length(td_use)
    subplot(3,4,4); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(1)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(1)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,8); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(2)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(2)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
    
    subplot(3,4,12); hold all;
    plot(td_use(t).(align_latent_signals_plot)(:,modesplot(3)),...
        'color',cols(t,:),'linewidth',3)
    plot(td_use(t).(align_latent_signals_plot)(1,modesplot(3)),...
        'color',cols(t,:),'marker','.','markersize',30)
    set(gca,'TickDir','out','FontSize',14), box off
end
set(f5,'color', [1 1 1]);



% set up unaligned Day 1
subplot(3,4,1);
axis tight;
set(gca,'YLim',xl);
ylabel(['Neural mode ' num2str(modesplot(1))]);
title(['Day 1'])
subplot(3,4,5);
axis tight;
set(gca,'YLim',yl);
ylabel(['Neural mode ' num2str(modesplot(2))]);
subplot(3,4,9);
axis tight;
set(gca,'YLim',zl);
ylabel(['Neural mode ' num2str(modesplot(3))]);


% set up unaligned Day N
subplot(3,4,2);
axis tight;
set(gca,'YLim',xl);
% ylabel(['Neural mode ' num2str(modesplot(1))]);
title(['Day ' num2str(diff_days)])
V = axis; text(0.5*V(2),0.9*V(3),['r=' num2str(r(1),3)],'FontSize',14);

subplot(3,4,6);
axis tight;
set(gca,'YLim',yl);
% ylabel(['Neural mode ' num2str(modesplot(2))]);
V = axis; text(0.5*V(2),0.9*V(3),['r=' num2str(r(2),3)],'FontSize',14);

subplot(3,4,10);
axis tight;
set(gca,'YLim',zl);
% ylabel(['Neural mode ' num2str(modesplot(3))]);
V = axis; text(0.5*V(2),0.9*V(3),['r=' num2str(r(3),3)],'FontSize',14);


% set up aligned Day 1
subplot(3,4,3);
axis tight;
set(gca,'YLim',xl);
% ylabel(['Neural mode ' num2str(modesplot(1))]);
title(['Day 1 Aligned'])

subplot(3,4,7);
axis tight;
set(gca,'YLim',yl);
% ylabel(['Neural mode ' num2str(modesplot(2))]);

subplot(3,4,11);
axis tight;
set(gca,'YLim',zl);
% ylabel(['Neural mode ' num2str(modesplot(3))]);



% set up aligned Day N
subplot(3,4,4);
axis tight;
set(gca,'YLim',xl);
% ylabel(['Neural mode ' num2str(modesplot(1))]);
title(['Day ' num2str(diff_days) ' Aligned'])
V = axis; text(0.5*V(2),0.9*V(3),['cc=' num2str(cc(1),3)],'FontSize',14);

subplot(3,4,8);
axis tight;
set(gca,'YLim',yl);
% ylabel(['Neural mode ' num2str(modesplot(2))]);
V = axis; text(0.5*V(2),0.9*V(3),['cc=' num2str(cc(2),3)],'FontSize',14);

subplot(3,4,12);
axis tight;
set(gca,'YLim',zl);
% ylabel(['Neural mode ' num2str(modesplot(3))]);
V = axis; text(0.5*V(2),0.9*V(3),['cc=' num2str(cc(3),3)],'FontSize',14);


%% bar graph of modes
f6 = figure('Position',[100 100 300 400]);
bar([r(1:3); cc(1:3)]');

h = legend({'Unaligned','Aligned'},'FontSize',14,'Location','South');
set(h,'Box','off');

set(gca,'Box','off','TickDir','out','FontSize',14,'XColor','k','YColor','k');
set(gca,'YLim',[-1 1]);
xlabel('Neural Mode');
ylabel('Correlation');

set(f6,'color', [1 1 1]);


%% visualize the weight matrices

cl = [min([min(min(A)) min(min(B))]), max([max(max(A)) max(max(B))])];

f7 = figure('Position',[100 100 800 400]);
subplot(1,2,1);
hold all;
imagesc(A);
set(gca,'YDir','reverse');
set(gca,'CLim',cl);
axis tight;
axis square;
set(gca,'Box','off','TickDir','out','FontSize',14);
title('A');
colorbar;

subplot(1,2,2);
hold all;
imagesc(B);
set(gca,'YDir','reverse');
set(gca,'CLim',cl);
axis tight;
axis square;
set(gca,'Box','off','TickDir','out','FontSize',14);
title('B');
colorbar;

set(f7,'color', [1 1 1]);



%%
if params.save_fig
    % Stability over time
    fn1 = [td(1).monkey '_' params.signals(1:end-4) '_Aligned_latent_trajectories_Day_1_' num2str(length(params.mani_dims)) 'D'];
    savefig(f1,fullfile(save_dir,params.signals(1:end-4),[fn1 '.fig']));
    saveas(f1,fullfile(save_dir,params.signals(1:end-4),[fn1 '.png']));
    saveas(f1,fullfile(save_dir,params.signals(1:end-4),[fn1 '.pdf']));
    
    fn2 = [td(1).monkey '_' params.signals(1:end-4) '_Aligned_latent_trajectories_Day_' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    savefig(f2,fullfile(save_dir,params.signals(1:end-4),[fn2 '.fig']));
    saveas(f2,fullfile(save_dir,params.signals(1:end-4),[fn2 '.png']));
    saveas(f2,fullfile(save_dir,params.signals(1:end-4),[fn2 '.pdf']));
    
    fn3 = [td(1).monkey '_' params.signals(1:end-4) '_Unaligned_latent_trajectories_Day_1_' num2str(length(params.mani_dims)) 'D'];
    savefig(f3,fullfile(save_dir,params.signals(1:end-4),[fn3 '.fig']));
    saveas(f3,fullfile(save_dir,params.signals(1:end-4),[fn3 '.png']));
    saveas(f3,fullfile(save_dir,params.signals(1:end-4),[fn3 '.pdf']));
    
    fn4 = [td(1).monkey '_' params.signals(1:end-4) '_Unaligned_latent_trajectories_Day_' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    savefig(f4,fullfile(save_dir,params.signals(1:end-4),[fn4 '.fig']));
    saveas(f4,fullfile(save_dir,params.signals(1:end-4),[fn4 '.png']));
    saveas(f4,fullfile(save_dir,params.signals(1:end-4),[fn4 '.pdf']));
    
    fn5 = [td(1).monkey '_' params.signals(1:end-4) '_Latent_trajectories_Against_Time_Day_1_to' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    savefig(f5,fullfile(save_dir,params.signals(1:end-4),[fn5 '.fig']));
    saveas(f5,fullfile(save_dir,params.signals(1:end-4),[fn5 '.png']));
    saveas(f5,fullfile(save_dir,params.signals(1:end-4),[fn5 '.pdf']));
    
    fn6 = [td(1).monkey '_' params.signals(1:end-4) '_Correlations_Day_1_to' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    savefig(f6,fullfile(save_dir,params.signals(1:end-4),[fn6 '.fig']));
    saveas(f6,fullfile(save_dir,params.signals(1:end-4),[fn6 '.png']));
    saveas(f6,fullfile(save_dir,params.signals(1:end-4),[fn6 '.pdf']));
    
    fn7 = [td(1).monkey '_' params.signals(1:end-4) '_Weight_Matrices_Day_1_to' num2str(diff_days) '_' num2str(length(params.mani_dims)) 'D'];
    savefig(f7,fullfile(save_dir,params.signals(1:end-4),[fn7 '.fig']));
    saveas(f7,fullfile(save_dir,params.signals(1:end-4),[fn7 '.png']));
    saveas(f7,fullfile(save_dir,params.signals(1:end-4),[fn7 '.pdf']));
    
end


