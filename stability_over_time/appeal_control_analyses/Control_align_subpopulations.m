%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  true;

monkey = 'chewie'; % 'chewie','chewie2','jaco','mihili'};
array = 'M1';
avg_dims = 1:4;
mani_dims = 1:10;


load([monkey '_data.mat']);

pars.mani_dims = mani_dims;
pars.align_latent_params.mani_dims = mani_dims;

dates = unique({master_td.date});
[align_results,decode_results,corr_results] = deal([]);
all_fr = [];
master_td_all_trials1 = [];
master_td_all_trials2 = [];

[r1, r2, tc1, tc2, good_idx] = deal([]);
for iDate = 1:length(dates)
    [~,td] = getTDidx(master_td,'date',dates{iDate});
    
    td = smoothSignals(td,struct('signals','M1_spikes','calc_fr',true));
    td2 =  td;
    
    idx =  randperm(size(td(1).M1_spikes,2));
    for trial = 1:length(td)
        td(trial).M1_spikes = td(trial).M1_spikes(:,idx(1:floor(length(idx)/2)));
        td2(trial).M1_spikes = td2(trial).M1_spikes(:,idx(floor(length(idx)/2)+1:end));
    end
    
    % compute PCA
    [td,~] = dimReduce(td,struct('signals','M1_spikes'));
    [td2,~] = dimReduce(td2,struct('signals','M1_spikes'));
    
    
    
    % compare original and nonlinear
    align_results =  [align_results, ...
        compDynamics( [td,td2], 'M1_pca', ...
        1:length(td), length(td)+1:length(td)+length(td2), pars.mani_dims )];
    
    corr_results = [corr_results, ...
        corrDynamics( [td,td2], 'M1_pca', ...
        1:length(td), length(td)+1:length(td)+length(td2), pars.mani_dims )];
    
    
    master_td_all_trials1 = [master_td_all_trials1, td];
    master_td_all_trials2 = [master_td_all_trials2, td2];

end

% compute cc metric for all days
[across,unalign] = deal(zeros(1,length(align_results)));
for iDate = 1:length(align_results)
    across(iDate) =  mean(align_results(iDate).cc(avg_dims));
    unalign(iDate) = mean(corr_results(iDate).r(avg_dims));
end


% do within day for real
within_align = align_latent_activity_within_day(master_td_all_trials, pars.align_latent_params );
[within1] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    within1(iDate) = mean(temp);
end


% do within day for sim
within_align = align_latent_activity_within_day(master_td_all_trials1, pars.align_latent_params );
[within2] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    within2(iDate) = mean(temp);
end



%% plot FR and tuning distributions

good_idx = good_idx == 1;
close all;

num_rows = 2;
num_cols = 2;

figure('Position',[100 100 1400 800]);


% plot example dynamics
subplot(num_rows,num_cols,1);
hold all;
plot(NaN,NaN,'b','LineWidth',3);
plot(NaN,NaN,'r','LineWidth',3);
trial = 1;
plot(td(trial).bin_size*(1:size(td(trial).M1_pca,1)), td(trial).M1_pca(:,1:4),'b','LineWidth',1)
plot(td2(trial).bin_size*(1:size(td2(trial).M1_pca,1)), td2(trial).M1_pca(:,1:4),'r','LineWidth',1);
set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('Latent Var. (a.u.)');
xlabel('Time (s)');
axis tight;
title('Top four latent variables')
h = legend({'Subpop 1','Subpop 2'},'Location','South');
set(h,'Box','off');



% plot CCs
subplot(num_rows,num_cols,2);
hold all;
% plot([1 1],[mean(real) - std(real), mean(real) + std(real)],'k-','LineWidth',2);
% plot(1,mean(real),'ko','LineWidth',2);
% plot([2 2],[mean(res) - std(res), mean(res) + std(res)],'k-','LineWidth',2);
% plot(2,mean(res),'ko','LineWidth',2)
boxplot([within2; within1; across; unalign]');

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'YLim',[0 1],'XLim',[0.5 4.5],'XTick',[1 2 3 4], ...
    'XTickLabel',{'Within Pop 1','Within Pop 2','Across','Unaligned'});
title('Canon. Corr.');
ylabel('Mean of top 4 CCs');


%

targs = unique([td.target_direction]);

plot_colors  =  distinguishable_colors(length(targs));

% get average for each trial
td_temp = td;
td2_temp = td2;
    
if 0
    td_temp = trimTD(td_temp,'idx_movement_on',{'idx_movement_on',5});
    td2_temp = trimTD(td2_temp,'idx_movement_on',{'idx_movement_on',5});
    td_temp = binTD(td_temp,'average');
    td2_temp = binTD(td2_temp,'average');
    
    ax(1) = subplot(num_rows,num_cols,3);
    hold all;
    for u = 1:length(targs)
        [~,temp] = getTDidx(td_temp,'target_direction',targs(u));
        temp = getSig(temp,{'M1_pca',1:3});
        plot3(temp(:,1),temp(:,2),temp(:,3),'.','MarkerSize',20,'Color',plot_colors(u,:));
    end
    axis square;
    
    set(gca,'Box','off','TickDir','out','FontSize',14);
    set(gca,'XLim',[-1.8, 1.8],'YLim',[-1.8 1.8],'ZLim',[-1.8 1.8]);
    set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
    title('Subpopulation 1');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    
    ax(2) = subplot(num_rows,num_cols,4);
    hold all;
    for u = 1:length(targs)
        [~,temp] = getTDidx(td2_temp,'target_direction',targs(u));
        temp = getSig(temp,{'M1_pca',1:3});
        plot3(temp(:,1),temp(:,2),temp(:,3),'.','MarkerSize',20,'Color',plot_colors(u,:));
    end
    axis square;
    
    set(gca,'Box','off','TickDir','out','FontSize',14);
    set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
    set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
    title('Subpopulation 2');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    set(ax(1),'CameraPosition',[8,10,15]);
    set(ax(2),'CameraPosition',[8,10,15]);
    linkprop(ax,'CameraPosition');
    
else
    td_temp = trialAverage(td_temp,'target_direction');
    td2_temp = trialAverage(td2_temp,'target_direction');
    
    ax(1) = subplot(num_rows,num_cols,3);
    hold all;
    for u = 1:length(targs)
        [~,temp] = getTDidx(td_temp,'target_direction',targs(u));
        temp = getSig(temp,{'M1_pca',1:3});
        plot3(temp(:,1),temp(:,2),temp(:,3),'-','LineWidth',2,'Color',plot_colors(u,:));
    end
    axis square;
    
    set(gca,'Box','off','TickDir','out','FontSize',14);
    set(gca,'XLim',[-1.8, 1.8],'YLim',[-1.8 1.8],'ZLim',[-1.8 1.8]);
    set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
    title('Subpopulation 1');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    
    ax(2) = subplot(num_rows,num_cols,4);
    hold all;
    for u = 1:length(targs)
        [~,temp] = getTDidx(td2_temp,'target_direction',targs(u));
        temp = getSig(temp,{'M1_pca',1:3});
        plot3(temp(:,1),temp(:,2),temp(:,3),'-','LineWidth',2,'Color',plot_colors(u,:));
    end
    axis square;
    
    set(gca,'Box','off','TickDir','out','FontSize',14);
    set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
    set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
    title('Subpopulation 2');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    set(ax(1),'CameraPosition',[8,10,15]);
    set(ax(2),'CameraPosition',[8,10,15]);
    linkprop(ax,'CameraPosition');
    
end


if save_figs
    saveas(gcf,fullfile(save_dir,'Subpopulation',[ monkey '_' array '_subpopulations.fig']));
    saveas(gcf,fullfile(save_dir,'Subpopulation',[ monkey '_' array '_subpopulations.pdf']));
    saveas(gcf,fullfile(save_dir,'Subpopulation',[ monkey '_' array '_subpopulations.png']));
end





