%% load data
clear all;
close all;
clc;


[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  false;


monkey = 'chewie'; % 'chewie','chewie2','jaco','mihili'};
avg_dims = 1:4;
pars.mani_dims = 1:10;
pars.align_latent_params.mani_dims = 1:10;


load([monkey '_data.mat']);

dates = unique({master_td.date});
[align_results,decode_results] = deal([]);
all_fr = [];
master_td_all_trials_sim = [];

[r1, r2, tc1, tc2, good_idx] = deal([]);
for iDate = 1:length(dates)
    [~,td] = getTDidx(master_td,'date',dates{iDate});
    
    td = smoothSignals(td,struct('signals','M1_spikes','calc_fr',true));
    
    % compute PCA
    [td,pca_info] = dimReduce(td,struct('signals','M1_spikes'));
    mu = mean(getSig(td,'M1_spikes'),1);
    
    % apply nonlinear transformation
    td2 = td;
    for trial = 1:length(td2)
        if 1
            td2(trial).M1_pca = td2(trial).M1_pca.* ...
                repmat(cos(linspace(-pi,pi,size((td2(trial).M1_pca),1)))', 1, size(td2(trial).M1_pca,2));
        else
            n = 10;
            for i = 1:size(td2(trial).M1_pca,2)
                td2(trial).M1_pca(:,i)  = conv(td2(trial).M1_pca(:,i),sin(linspace(-pi,pi,n)),'same');
            end
        end
        
        % project back into firing rate space
        td2(trial).M1_spikes = td2(trial).M1_pca * inv(pca_info.w) + repmat(mu,size(td2(trial).M1_spikes,1),1);
    end
    
    all_fr = [all_fr, [mean(getSig(td,'M1_spikes'),1); mean(getSig(td2,'M1_spikes'),1)] ];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fit tuning curves
    fr = getSig(td,'M1_spikes');
    v = getSig(td,'vel');
    s = sqrt(v(:,1).^2 + v(:,2).^2);
    
    th = atan2(v(:,2),v(:,1));
    X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
    
    for unit = 1:size(fr,2)
        [b,~,~,~,temp] = regress(fr(:,unit),X);
        r1 = [r1; temp(1)];
        % convert to model b0 + b1*cos(theta+b2)
        tc1  = [tc1; b(1), sqrt(b(2).^2 + b(3).^2), atan2(b(2),b(3))];
        good_idx = [good_idx; temp(1) > 0.3];
    end
    
    
    fr = getSig(td2,'M1_spikes');
    v = getSig(td2,'vel');
    s = sqrt(v(:,1).^2 + v(:,2).^2);
    
    th = atan2(v(:,2),v(:,1));
    X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
    
    for unit = 1:size(fr,2)
        [b,~,~,~,temp] = regress(fr(:,unit),X);
        r2 = [r2; temp(1)];
        % convert to model b0 + b1*cos(theta+b2)
        tc2  = [tc2; b(1), sqrt(b(2).^2 + b(3).^2), atan2(b(2),b(3))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % compare original and nonlinear
    align_results =  [align_results, ...
        compDynamics( [td,td2], 'M1_pca', ...
        1:length(td), length(td)+1:length(td)+length(td2), pars.mani_dims )];
    
    
    % do decoding
    pars.decoder_params.hist_bins = 1;
    for trial = 1:length(td2)
        td2(trial).date = datestr(datenum(td2(trial).date,'mm-dd-yyyy')+1,'mm-dd-yyyy');
    end
    decode_results = [decode_results, ...
        decode_across_days([td td2],pars.decoder_params)];
    
    
    
    
    % now build master vector with all trials with the sim
    [~,td] = getTDidx(master_td_all_trials,'date',dates{iDate});
    td = smoothSignals(td,struct('signals','M1_spikes','calc_fr',true));
    % compute PCA
    [td,pca_info] = dimReduce(td,struct('signals','M1_spikes'));
    mu = mean(getSig(td,'M1_spikes'),1);
    % apply nonlinear transformation
    td2 = td;
    for trial = 1:length(td2)
        td2(trial).M1_pca = td2(trial).M1_pca.*repmat(cos(linspace(-pi,pi,size((td2(trial).M1_pca),1)))',1,size(td2(trial).M1_pca,2));
        
        % project back into firing rate space
        td2(trial).M1_spikes = td2(trial).M1_pca * inv(pca_info.w) + repmat(mu,size(td2(trial).M1_spikes,1),1);
    end
    
    td2 = dimReduce(td2,'M1_spikes');
    
    master_td_all_trials_sim = [master_td_all_trials_sim, td2];
end

% compute cc metric for all days
[across] = deal(zeros(1,length(align_results)));
for iDate = 1:length(align_results)
    across(iDate) =  mean(align_results(iDate).cc(avg_dims));
end


% do within day for real
within_align = align_latent_activity_within_day(master_td_all_trials, pars.align_latent_params );
[within_real] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    within_real(iDate) = mean(temp);
end


% do within day for sim
within_align = align_latent_activity_within_day(master_td_all_trials_sim, pars.align_latent_params );
[within_sim] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    within_sim(iDate) = mean(temp);
end



% compute decoding results
dec_within = cellfun(@mean,{decode_results.withinR2_m});
dec_across = cellfun(@mean,{decode_results.acrossR2});
dec_ctrl = cellfun(@mean,{decode_results.ctrlR2});


%% plot FR and tuning distributions

good_idx = good_idx == 1;
close all;

num_rows = 2;
num_cols = 4;

figure('Position',[100 100 1400 800]);


subplot(num_rows,num_cols,6);
hold all;
boxplot(all_fr'./td(1).bin_size);
set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('Firing Rate (Hz)');
title('Mean Single Neuron FR');
set(gca,'XTick',[1 2],'XTickLabel',{'Real Day 1','Sim Day 2'});



% plot PDs
subplot(num_rows,num_cols,7);
hold all;
[n,x] = hist(tc1(:,3),[-pi:pi/8:pi]);
h = bar(x,n,'hist');
set(h,'FaceColor','b','FaceAlpha',0.3);
[n,x] = hist(tc2(:,3),[-pi:pi/8:pi]);
h = bar(x,n,'hist');
set(h,'FaceColor','r','FaceAlpha',0.3);
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XLim',[-pi pi]);
title('Distribution of Single Neuron PDs');
xlabel('PD (rad)');
ylabel('Count');
h = legend({'Real Day 1','Sim Day 2'});
set(h,'Box','off');

subplot(num_rows,num_cols,8);
hold all;
[n,x] = hist(angleDiff(tc1(:,3),tc2(:,3),true,true),[-pi:pi/8:pi]);
h = bar(x,n,'hist');
set(h,'FaceColor','k','FaceAlpha',0.3);
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XLim',[-pi pi]);
title('Diff between Real and Sim Neurons');
xlabel('Change in PD (rad)');
ylabel('Count');



% plot example dynamics
subplot(num_rows,num_cols,2);
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
h = legend({'Real Day 1','Sim Day 1'},'Location','South');
set(h,'Box','off');



% plot CCs
subplot(num_rows,num_cols,3);
hold all;
% plot([1 1],[mean(real) - std(real), mean(real) + std(real)],'k-','LineWidth',2);
% plot(1,mean(real),'ko','LineWidth',2);
% plot([2 2],[mean(res) - std(res), mean(res) + std(res)],'k-','LineWidth',2);
% plot(2,mean(res),'ko','LineWidth',2)
boxplot([within_sim; within_real; across]');

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'YLim',[0 1],'XLim',[0.5 3.5],'XTick',[1 2 3], ...
    'XTickLabel',{'Within Sim','Within Real','Across (Sim vs Real)'});
title('Canon. Corr.');
ylabel('Mean of top 4 CCs');



% plot decoding
subplot(num_rows,num_cols,4);
hold all;
boxplot([dec_ctrl; dec_within; dec_across]');
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',1:3, ...
    'XTickLabel',{'Within Sim','Within  Real','Across (Sim vs Real)'});
set(gca,'YLim',[0 1]);
ylabel('R^2');
title('Decoder Perf.');


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
    
    ax(1) = subplot(2,4,1);
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
    title('Real Day 1 Population');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    
    ax(2) = subplot(2,4,5);
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
    title('Simulated Day 2 Population');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    set(ax(1),'CameraPosition',[8,10,15]);
    set(ax(2),'CameraPosition',[8,10,15]);
    linkprop(ax,'CameraPosition');
else
    td_temp = trialAverage(td_temp,'target_direction');
    td2_temp = trialAverage(td2_temp,'target_direction');
    
    ax(1) = subplot(2,4,1);
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
    title('Real Day 1 Population');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    
    ax(2) = subplot(2,4,5);
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
    title('Simulated Day 2 Population');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    set(ax(1),'CameraPosition',[8,10,15]);
    set(ax(2),'CameraPosition',[8,10,15]);
    linkprop(ax,'CameraPosition');
    
end


if save_figs
    saveas(gcf,fullfile(save_dir,'NonlinearSim',[ monkey '_nonlinear-sim.fig']));
    saveas(gcf,fullfile(save_dir,'NonlinearSim',[ monkey '_nonlinear-sim.pdf']));
    saveas(gcf,fullfile(save_dir,'NonlinearSim',[ monkey '_nonlinear-sim.png']));
end




