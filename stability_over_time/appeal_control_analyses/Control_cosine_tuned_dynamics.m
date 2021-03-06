% divide by 1/sqrt(n)
% population with random 0,-1,+1 such taht the sum across modes is the same
% for all
% we want the modes to be orthogonal from neural space
% the sum across modes of weights has to be the same




%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  false;

do_speed = false;

monkey = 'chewie'; % 'chewie','chewie2','jaco','mihili'};
array = 'M1';

avg_dims = 1:4;
pars.mani_dims = 1:8;
pars.align_latent_params.mani_dims = 1:8;


load([monkey '_data.mat']);

dates = unique({master_td.date});
[align_results,decode_results,decode_results_unaligned,corr_results] = deal([]);
all_fr = [];
master_td_all_trials_sim = [];

v_max = max(prctile(getSig(getNorm(master_td,'vel'),'vel_norm'),[2.5 99.9]));

[r1, r2, tc1, tc2, good_idx] = deal([]);
for iDate = 1:length(dates)
    [~,td] = getTDidx(master_td,'date',dates{iDate});
    
    td = smoothSignals(td,struct('signals','M1_spikes','calc_fr',true));
    mu = mean(getSig(td,'M1_spikes'));
    
    targs = unique([td.target_direction]);
    
    % generate tuning curves for the modes
    n = length(targs);
    tc = [zeros(n,1), sin(targs)', cos(targs)', 0.1*ones(n,1)./v_max];
    
    
    % generate dynamics that are cosine tuned w.r.t. kinematics
    scale_fac = 0.05;
    % sort by "variance explained"
    [~,I] = sort(abs(tc(:,1)+tc(:,2)+tc(:,3)));
    tc = tc(I,:);
    
    
    % generate fake PC matrix
    vals = { ...
        [1 1 1 1 1 1 1 1], ...
        [1 0 1 0 1 0 1 0], ...
        [0 1 0 1 0 1 0 1], ...
        [1 0 0 0 1 0 0 0], ...
        [0 1 0 0 0 1 0 0], ...
        [0 0 1 0 0 0 1 0], ...
        [0 0 0 1 0 0 0 1], ...
        };
    
    w = zeros(n,size(td(1).M1_spikes,2));
    for i = 1:size(td(1).M1_spikes,2)
        w(:,i) = (-1+2*rand)*vals{randi(7)}';
        %w(i,:) = vals(randperm(n));
    end
    
    w = w./sqrt(n);
    
    % loop  along directions
    td2 = td;
    for trial = 1:length(td)
        v = getSig(td(trial),'vel');
        s = sqrt(v(:,1).^2 + v(:,2).^2);
        
        th = atan2(v(:,2),v(:,1));
        
        lat_vars = zeros(size(s,1),n);
        for i = 1:n
            lat_vars(:,i) = scale_fac * abs(tc(i,1) + tc(i,4)*s + tc(i,2).*s.*sin(th) + tc(i,3).*s.*cos(th));
        end
        
        temp = lat_vars * w + repmat(mu,size(td2(trial).M1_spikes,1),1) + normrnd(0,0.005,size(td2(trial).M1_spikes));
        temp(temp < 0) = 0;
        %         temp(temp > 10) = 10;
        td2(trial).M1_spikes = temp;
    end
    
    td2 = smoothSignals(td2,struct('signals','M1_spikes','calc_fr',true));
    
    
    % compute PCA
    [td,pca_info] = dimReduce(td,struct('signals','M1_spikes'));
    [td2,~] = dimReduce(td2,struct('signals','M1_spikes'));
    
    all_fr = [all_fr, [mean(getSig(td,'M1_spikes'),1); mean(getSig(td2,'M1_spikes'),1)] ];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fit tuning curves
    fr = getSig(td,'M1_spikes');
    v = getSig(td,'vel');
    s = sqrt(v(:,1).^2 + v(:,2).^2);
    
    th = atan2(v(:,2),v(:,1));
    
    if do_speed
        X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
    else
        X = [ones(size(v,1),1) sin(th) cos(th)];
    end
    
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
    
    if do_speed
        X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
    else
        X = [ones(size(v,1),1) sin(th) cos(th)];
    end
    
    
%     figure; hold all;
%     for unit = 1:size(fr,2)
%         
%         plot3(X(:,2),X(:,3),fr(:,unit),'.')
%         pause;
%         %         close all;
%     end
    
    
    
    for unit = 1:size(fr,2)
        [b,~,~,~,temp] = regress(fr(:,unit),X);
        r2 = [r2; temp(1)];
        % convert to model b0 + b1*cos(theta+b2)
        tc2  = [tc2; b(1), sqrt(b(2).^2 + b(3).^2), atan2(b(2),b(3))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % compare original and shuffled
    align_results =  [align_results, ...
        compDynamics( [td,td2], 'M1_pca', ...
        1:length(td), length(td)+1:length(td)+length(td2), pars.mani_dims )];
    
    % compare unaligned
    corr_results = [corr_results, ...
        corrDynamics( [td,td2], 'M1_pca', ...
        1:length(td), length(td)+1:length(td)+length(td2), pars.mani_dims )];
    
    
    % do decoding
    pars.decoder_params.hist_bins = 1;
    for trial = 1:length(td2)
        td2(trial).date = datestr(datenum(td2(trial).date,'mm-dd-yyyy')-1,'mm-dd-yyyy');
    end
    pars.decoder_params.in =  'aligned_data';
    decode_results = [decode_results, ...
        decode_across_days([td td2],pars.decoder_params)];
    
    pars.decoder_params.in =  'unaligned_data';
    decode_results_unaligned = [decode_results_unaligned, ...
        decode_across_days([td td2],pars.decoder_params)];
    
    
    master_td_all_trials_sim = [master_td_all_trials_sim, td2];
end

% compute cc metric for all days
[across,corr] = deal(zeros(1,length(align_results)));
for iDate = 1:length(align_results)
    across(iDate) =  mean(align_results(iDate).cc(avg_dims));
    corr(iDate) = mean(corr_results(iDate).r(avg_dims));
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
dec_unaligned = cellfun(@mean,{decode_results_unaligned.acrossR2});


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
ax(1) = subplot(num_rows,num_cols,7);
hold all;
[n,x] = hist(r1,0:0.1:1);
h = bar(x,n,'hist');
set(h,'FaceColor','k','FaceAlpha',0.3);
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XLim',[0 1]);
title('Real Neurons');
xlabel('R^2');
ylabel('Count');

ax(2) = subplot(num_rows,num_cols,8);
hold all;
[n,x] = hist(r2,0:0.1:1);
h = bar(x,n,'hist');
set(h,'FaceColor','k','FaceAlpha',0.3);
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XLim',[0 1]);
title('Simulation');
xlabel('R^2');
ylabel('Count');

V1 = axis(ax(1));
V2 = axis(ax(2));
set(ax(1),'YLim',[0  max([V1(4) V2(4)])]);
set(ax(2),'YLim',[0  max([V1(4) V2(4)])]);


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
boxplot([within_sim; across; corr]');

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'YLim',[0 1],'XLim',[0.5 3.5],'XTick',[1 2 3], ...
    'XTickLabel',{'Within Shuffle','Across (Aligned)','Across (Unaligned)'}, ...
    'XTickLabelRotation',30);
title('Canon. Corr.');
ylabel('Mean of top 4 CCs');



% plot decoding
subplot(num_rows,num_cols,4);
hold all;
boxplot([dec_within; dec_across; dec_unaligned]');
set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',1:3, ...
    'XTickLabel',{'Within','Across (Aligned)','Across (Unaligned)'}, ...
    'XTickLabelRotation',30);
set(gca,'YLim',[0 1]);
ylabel('R^2');
title('Decoder Perf.');


%

targs = unique([td.target_direction]);

plot_colors  =  distinguishable_colors(length(targs));

% get average for each trial
td_temp = td;
td2_temp = td2;

if 1
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
    %     set(gca,'XLim',[-1.8, 1.8],'YLim',[-1.8 1.8],'ZLim',[-1.8 1.8]);
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
    %     set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
    set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
    title('Simulated Day 2 Population');
    
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    
    %     set(ax(1),'CameraPosition',[8,10,15]);
    %     set(ax(2),'CameraPosition',[8,10,15]);
    %     linkprop(ax,'CameraPosition');
    
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
    saveas(gcf,fullfile(save_dir,'Cosine Tuned Dynamics',[ monkey '_' array '_CosineTunedDynamics.fig']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuned Dynamics',[ monkey '_' array '_CosineTunedDynamics.pdf']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuned Dynamics',[ monkey '_' array '_CosineTunedDynamics.png']));
end



