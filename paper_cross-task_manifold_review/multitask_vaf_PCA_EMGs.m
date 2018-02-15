%
% Multi-task VAF EMG modes
%


sessions_to_use = [7 8 9 1:3]; % wrist: [7 8 9 1:3] % reach: [4:6 10:11]


% preallocate PCA VAF matrix
max_n_emgs = max( cellfun(@(x) length(x.chosen_emgs), datasets(sessions_to_use)) );
vaf_emg = nan(length(sessions_to_use),max_n_emgs); 


for s = 1:length(sessions_to_use)
    
    n_tasks = length(datasets{sessions_to_use(s)}.labels);    

    stdata = datasets{sessions_to_use(s)}.stdata;
    
    % get concatenated single trials for all tasks in this session
    multitask_emg = cell2mat(cellfun(@(x) x.target{end}.emg_data.conc_emg, stdata,'UniformOutput',false)');
    
    % Do PCA
    multitask_vaf = svd(cov(multitask_emg))/sum(svd(cov(multitask_emg)));
    
    vaf_emg(s,1:length(multitask_vaf)) = multitask_vaf;    
end



n_dims_95vaf = zeros(1,length(sessions_to_use));

% Find number of dimensions that explain >= 95 % of the variance
for s = 1:length(sessions_to_use)
    n_dims_95vaf(s) = find(cumsum(vaf_emg(s,:),2)>.95,1);
end


% Plot VAF for all tasks from one session --raw figure
cols = parula(size(vaf_emg,1)+1);

figure,hold on
for s = 1:length(sessions_to_use)
    t_data = cumsum(vaf_emg(s,~isnan(vaf_emg(s,:))))*100;
    plot(t_data,'color',cols(s,:))
    plot(n_dims_95vaf(s),t_data(n_dims_95vaf(s)),'.','markersize',32,'color',cols(s,:))
    xlabel('EMG modes'), xlim([0 max_n_emgs])
    set(gca,'TickDir','out','FontSize',14), box off
    ylim([0 100])
    if s == 1, ylabel('EMG variance expl. (%)'); end
end


% plot histogram
x_ax = 1:max_n_emgs+1;
y_hist = histcounts(n_dims_95vaf,x_ax);
mn_hist = mean(n_dims_95vaf);
sd_hist = std(n_dims_95vaf);

figure,hold on
bar(x_ax(1:end-1),y_hist,'FaceColor',[.6 .6 .6],'Edgecolor',[.6 .6 .6])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('EMG modes >= 95 % multi-task variance'), xlim([0 max_n_emgs]), ylabel('Counts')
yl = ylim; 
plot(mn_hist,yl(2)+.3,'.','markersize',40,'color',[.6 .6 .6])
plot([mn_hist-sd_hist, mn_hist+sd_hist],[yl(2)+.3,yl(2)+.3],'linewidth',1.5,'color',[.6 .6 .6])