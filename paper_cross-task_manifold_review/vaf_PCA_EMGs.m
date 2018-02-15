%
% Compute VAF of the EMGs using PCA. This figure generates Figure 2 in the
% Nature Comms review letter
%


sessions_to_use = [7 8 9 1:3]; % wrist: [7 8 9 1:3] % reach: [4:6 10:11]


for s = 1:length(sessions_to_use)
    
    d_pca_emg = sessions_to_use(s);

    n_tasks = length(datasets{d_pca_emg}.labels);
    vaf_emg = zeros(n_tasks,length(datasets{d_pca_emg}.chosen_emgs)); % preallocate PCA VAF matrix

    % Do PCA on the trial-related part of the data
    for t = 1:n_tasks

        emg = datasets{d_pca_emg}.stdata{t}.target{end}.emg_data.conc_emg;    
        vaf_emg(t,:) = svd(cov(emg))/sum(svd(cov(emg)));
    end
    
    
%     % Plot VAF for all tasks from one session --raw figure
%     cols = [.4 .6 .5; .7 .1 .1; .5 .65 .9; 1 .6 .3];
%     
%     figure,
%     for t = 1:n_tasks
%         subplot(1,n_tasks,t), hold on
%         bar(cumsum(vaf_emg(t,:))*100,'FaceColor',cols(t,:),'EdgeColor',cols(t,:))
%         plot([0 size(emg,2)+1],[95 95],'linewidth',2,'color',[.6 .6 .6])
%         xlabel('EMG modes'), xlim([0 size(emg,2)+1])
%         set(gca,'TickDir','out','FontSize',14), box off
%         if t == 1, ylabel('EMG variance expl. (%)'); end
%         title(datasets{d_pca_emg}.labels{t},'Fontsize',14)
%     end


    % save in struct
    pca_emg{s} = vaf_emg;
end



% -------------------------------------------------------------------------
% Plot VAF for all tasks from one session
ex_session = 2;
n_tasks = length(datasets{sessions_to_use(ex_session)}.labels);
n_emgs = length(datasets{sessions_to_use(ex_session)}.chosen_emgs);
cols = [.4 .6 .5; .7 .1 .1; .5 .65 .9; 1 .6 .3];

figure,
for t = 1:n_tasks
    subplot(1,n_tasks,t), hold on
    bar(cumsum(pca_emg{ex_session}(t,:))*100,'FaceColor',cols(t,:),'EdgeColor',cols(t,:))
    plot([0 n_emgs+1],[95 95],'linewidth',2,'color',[.6 .6 .6])
    xlabel('EMG modes'), xlim([0 n_emgs+1])
    set(gca,'TickDir','out','FontSize',14), box off
    if t == 1, ylabel('EMG variance expl. (%)'); end
    title(datasets{sessions_to_use(ex_session)}.labels{t},'Fontsize',14)
end



% -------------------------------------------------------------------------
% Histogram of number of dimensions to explain more than some percentage of
% the EMG variance
perc_var = 95;
x_hist = 1:10;

n_dims_perc_var = [];

for s = 1:length(sessions_to_use)
    for t = 1:size(pca_emg{s},1)
       n_dims_perc_var = [n_dims_perc_var, find(cumsum(pca_emg{s}(t,:))>perc_var/100,1)]; 
    end
end

hist_n_dims = histcounts(n_dims_perc_var,x_hist);
yhistplot = max(hist_n_dims/sum(hist_n_dims)*100);
figure,bar(x_hist(1:end-1),hist_n_dims/sum(hist_n_dims)*100,'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off
xlabel(['Nbr. EMG modes explain >' num2str(perc_var) ' % variance'])
ylabel('Percentage (%)')
text(x_hist(end)-2,yhistplot+10,['n=' num2str(sum(n_dims_perc_var))],'FontSize',16)
ylim([0 yhistplot+20])
hold on, plot(mean(n_dims_perc_var),yhistplot+10,'.','markersize',24,'color',[.6 .6 .6])
plot([mean(n_dims_perc_var)-std(n_dims_perc_var), mean(n_dims_perc_var)+std(n_dims_perc_var)],...
    [yhistplot+10, yhistplot+10],'linewidth',3,'color',[.6 .6 .6])

% % Plot VAF for all the tasks on top of each other
% figure, hold on
% for s = 1:length(sessions_to_use)
%    
%     for t = 1:size(pca_emg{s},1) % color code by monkey
%         if s<=3
%             plot(cumsum(pca_emg{s}(t,:)),'color','k');
%         else
%             plot(cumsum(pca_emg{s}(t,:)),'color',[.6 .6 .6]);
%         end
%     end
%     xlabel('EMG modes'), ylabel('EMG variance expl. (%)'),ylim([0 1])
%     set(gca,'TickDir','out','FontSize',14), box off
% end



clearvars -except datasets
