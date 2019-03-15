%
% TME control analyses to check whether the stable population dynamics
% arise from features in the single unit data. Many options!
%
%



rng('shuffle', 'twister') % randomize the seed


% TME parameters
surrogate_type = 'surrogate-TdiagC';% 'surrogate-T+-1'; % 'surrogate-TdiagC'; % 'surrogate-NC'; 'surrogate_T'; 'surrogate_TN'; 'surrogate_C'; 
n_surrogates = 2000; % increase to 10000

plot_all_yn = false;




% do for all sessions
for s = 1:n_sessions
   
    % get trials for this session
    [~, td] = getTDidx(master_td,{'date',meta.sessions{s}});

    
    % trial average the data, as required by TME
    tda = trialAverage(td,{'target_direction'});
    
    
    % PUT THE DATA INTO THE RIGHT FORMAT: FR: time x unit x condition --I think
    n_targets = length(unique(meta.targets));
    
    fr = zeros( size(td(1).(pars.spiking_inputs{1}),1), size(td(1).(pars.spiking_inputs{1}),2), n_targets);
    
    % fill the FR tensor
    for t = 1:n_targets   
        fr(:,:,t) = tda(t).(pars.spiking_inputs{1});
    end

    
    % -------------------------------------------------------------------------
    %% GET PRIMARY FEATURES OF ORIGINAL DATA
    
    % Quantify primary features of the original data
    [targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures_J(fr);
    
    
    % define surrogate params (T,N,C)
    if strcmp(surrogate_type, 'surrogate-T')
        params.margCov{1} = targetSigmaT;
        params.margCov{2} = [];
        params.margCov{3} = [];
        params.meanTensor = M.T;
    elseif strcmp(surrogate_type, 'surrogate-TdiagC')
        params.margCov{1} = targetSigmaT;
        params.margCov{2} = [];
        params.margCov{3} = eye(size(targetSigmaC)).*diag(targetSigmaC);
        params.meanTensor = M.TC;
    elseif strcmp(surrogate_type, 'surrogate-TN')
        params.margCov{1} = targetSigmaT;
        params.margCov{2} = targetSigmaN;
        params.margCov{3} = [];
        params.meanTensor = M.TN;
    elseif strcmp(surrogate_type, 'surrogate-T+-1')
        targetSigmaTpdiag = targetSigmaT;
        for i = 1:size(targetSigmaTpdiag,1)
            targetSigmaTpdiag(i,i+2:end) = 0;
            if i > 2
                targetSigmaTpdiag(i,1:i-2) = 0;
            end
        end
        params.margCov{1} = targetSigmaTpdiag;
        params.margCov{2} = [];
        params.margCov{3} = [];
        params.meanTensor = M.T;
    elseif strcmp(surrogate_type, 'surrogate-NC')
        params.margCov{1} = [];
        params.margCov{2} = targetSigmaN;
        params.margCov{3} = targetSigmaC;
        params.meanTensor = M.NC;
    elseif strcmp(surrogate_type, 'surrogate-C')
        params.margCov{1} = [];
        params.margCov{2} = [];
        params.margCov{3} = targetSigmaC;
        params.meanTensor = M.C;
    elseif strcmp(surrogate_type, 'surrogate-TNC')
        params.margCov{1} = targetSigmaT;
        params.margCov{2} = targetSigmaN;
        params.margCov{3} = targetSigmaC;
        params.meanTensor = M.TNC;
    else
        error('Need to code these TME parameters')
    end
    
    
    % Fit the maximum entropy distribution
    maxEntropy = fitMaxEntropy(params);
    
    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    %% COMPUTE PRINCIPAL ANGLES BETWEEN SURROGATE AND ACTUAL DATASETS
    
    surrCCs = zeros(n_surrogates,length(pars.mani_dims));
    pca_signals = [pars.spiking_inputs{1}(1:end-7) '_pca'];
    
    for su = 1:n_surrogates
        
        if rem(su,100) == 0
            fprintf('surrogate %d from %d \n', su, n_surrogates)
        end
        
        % 1. generate surrogate datasets
        surr_tensor = sampleTME(maxEntropy);
        
        % get covariances
        [surrSigmaT, surrSigmaN, surrSigmaC, Msurr] = extractFeatures_J(surr_tensor);

        
        % 2. Do PCA OF THE SURROGATE DATASETS
        surr_fr = permute(surr_tensor,[2 1 3]);
        
        surr_fr = reshape(surr_fr,size(surr_fr,1),[]);
        
        [lvsu, wsu] = pca(surr_fr);
        
        % 3. Retrieve actual data into a matrix with the same size
        lva = [];
        for t = 1:n_targets
           
            lva = cat(1,lva,td(t).(pca_signals));
        end
        
        % 4. Keep the manifold dimensions we care about
        lva = lva(:,pars.mani_dims);
        lvsu = lvsu(:,pars.mani_dims);
        
        
        % 5. do CCA to align the dynamics
        [~, ~, cc] = canoncorr(lva,lvsu);
        
        surrCCs(su,:) = cc;
    end

    
    % SAVE DATA FOR THIS SESSION
    
    
    % Retrieve within-day data
    wdCCs = zeros(length(within_day_align_results(s).aligned_info),length(pars.mani_dims));
    for w = 1:length(within_day_align_results(s).aligned_info)
        wdCCs(w,:) = within_day_align_results(s).aligned_info(w).cc;
    end
    
    TME_control(s).surrCCs = surrCCs';
    TME_control(s).wsCCs = wdCCs';
    


    % PLOT PER DATASET
    if plot_all_yn
    
        % GET SUMMARY STATISTICS
        
        mn_cc_su = mean(surrCCs(:,1:pars.align_latent_params.top_lv_plot),2);
        
        x_hist = 0:0.05:1.05;
        y_cc_su = histcounts(mn_cc_su,x_hist)/numel(mn_cc_su)*100;
        
        mn_cc_wd = mean(wdCCs(:,1:pars.align_latent_params.top_lv_plot),2);
        y_cc_wd = histcounts(mn_cc_wd,x_hist)/numel(mn_cc_wd)*100;
        
        
        
        % PLOT THAT COMPARES TO WITHIN-DAY
        figure,hold on
        hsu = bar(x_hist(1:end-1),y_cc_su,'hist');
        hswd = bar(x_hist(1:end-1),y_cc_wd,'hist');
        set(hsu, 'FaceColor',[.6 .6 .6])
        alpha(hsu,.5)
        alpha(hswd,.5)
        set(hswd, 'FaceColor','r')
        xlim([0 1]);
        set(gcf,'color','w')
        set(gca,'TickDir','out','FontSize',14)
        legend('surrogate same tuning','actual','Location','NorthWest'), legend boxoff
        xlabel(['Mean CC top ' num2str(pars.align_latent_params.top_lv_plot) ' modes']),ylabel('Probability (%)')
    end
end




% POOL DATA ACROS SESSSIONS FOR SUMMARY HISTOGRAM
pooled_surrCCs = [TME_control.surrCCs]';
pooled_withinCCs = [TME_control.wsCCs]';




%  ========================================================================
%% SUMMARY PLOT
mn_cc_su = mean(pooled_surrCCs(:,1:pars.align_latent_params.top_lv_plot),2);

x_hist = 0:0.05:1.05;
y_cc_su = histcounts(mn_cc_su,x_hist)/numel(mn_cc_su)*100;

mn_cc_wd = mean(pooled_withinCCs(:,1:pars.align_latent_params.top_lv_plot),2);
y_cc_wd = histcounts(mn_cc_wd,x_hist)/numel(mn_cc_wd)*100;


mn_mn_cc_wd = mean(mn_cc_wd);
sd_mn_cc_wd = std(mn_cc_wd);
mn_mn_cc_su = mean(mn_cc_su);
sd_mn_cc_su = std(mn_cc_su);


% PLOT THAT COMPARES TO WITHIN-DAY
figure,hold on
hsu = bar(x_hist(1:end-1),y_cc_su,'hist');
hswd = bar(x_hist(1:end-1),y_cc_wd,'hist');
set(hsu, 'FaceColor',[.6 .6 .6])
alpha(hsu,.5)
alpha(hswd,.5)
set(hswd, 'FaceColor','r')
xlim([0 1]);
set(gcf,'color','w')
set(gca,'TickDir','out','FontSize',14)
xlabel(['Mean CC top ' num2str(pars.align_latent_params.top_lv_plot) ' modes']),ylabel('Probability (%)')
title(['Pooled data - ' pars.monkey ' - ' pars.spiking_inputs{1}(1:end-7)])
yl = ylim;
ylim([0 yl(2)+yl(2)/5])
yax = yl(2)+yl(2)/10;
plot(mn_mn_cc_wd,yax,'.','color','r','markersize',24)
plot([mn_mn_cc_wd-sd_mn_cc_wd,mn_mn_cc_wd+sd_mn_cc_wd],[yax yax],'color','r','linewidth',1.5)
plot(mn_mn_cc_su,yax,'.','color',[.6 .6 .6],'markersize',24)
plot([mn_mn_cc_su-sd_mn_cc_su,mn_mn_cc_su+sd_mn_cc_su],[yax yax],'color',[.6 .6 .6],'linewidth',1.5)
legend(surrogate_type,'actual','Location','NorthWest'), legend boxoff

%  ========================================================================
%% PLOT SOME RAW DATA

cols = parula(n_targets+1);


% PLOT RANDOM SURROGATE NEURONS

% rnd_surr_units = randperm(size(surr_tensor,2));
% rnd_surr_units = rnd_surr_units(1:2);
% 
% figure
% subplot(221), hold on
% for t = 1:n_targets
%     plot(surr_tensor(:,rnd_surr_units(1),t),'linewidth',1,'color',cols(t,:));
% end
% subplot(222), hold on
% for t = 1:n_targets
%     plot(surr_tensor(:,rnd_surr_units(2),t),'linewidth',1,'color',cols(t,:));
% end



% PLOT CLOUD POINTS!

% put the actual data into tensor form (time x target x LV)
actual_tensor = reshape(lva,size(lva,1)/n_targets,n_targets,size(lva,2));

modes_plot = 2:4;
bin_start = 8;
bin_end = size(actual_tensor,1);

figure,
subplot(121), hold on
for t = 1:n_targets
    plot3(squeeze(actual_tensor(bin_start:bin_end,t,modes_plot(1))),...
        squeeze(actual_tensor(bin_start:bin_end,t,modes_plot(2))),...
        squeeze(actual_tensor(bin_start:bin_end,t,modes_plot(3))),...
        'color',cols(t,:),'marker','.','markersize',16,'linestyle','none');
end
title('Actual data')

subplot(122), hold on
for t = 1:n_targets
    plot3(squeeze(surr_tensor(bin_start:bin_end,t,modes_plot(1))),...
        squeeze(surr_tensor(bin_start:bin_end,t,modes_plot(2))),...
        squeeze(surr_tensor(bin_start:bin_end,t,modes_plot(3))),...
        'color',cols(t,:),'marker','.','markersize',16,'linestyle','none');
end
title('Surrogate data')


% PLOT COVARIANCE MATRICES ACTUAL AND SURROGATE DATA, and scatter plot
% comparison beween the two
% 
% figure('units','normalized','outerposition',[0.2 0.1 .6 .8])
% subplot(331),
% imagesc(targetSigmaT)
% set(gca,'TickDir','out','FontSize',14)
% title('Cov Time - Actual')
% axis square
% 
% subplot(334),
% imagesc(targetSigmaN)
% set(gca,'TickDir','out','FontSize',14)
% title('Cov Neurons - Actual')
% axis square
% 
% subplot(337),
% imagesc(targetSigmaC)
% set(gca,'TickDir','out','FontSize',14)
% title('Cov Cond. - Actual')
% axis square
% 
% subplot(332),
% imagesc(surrSigmaT)
% set(gca,'TickDir','out','FontSize',14)
% ylabel('Cov Time')
% axis square
% title('Cov Time - Surr.')
% 
% subplot(335),
% imagesc(surrSigmaN)
% set(gca,'TickDir','out','FontSize',14)
% title('Cov Neurons - Surr.')
% axis square
% 
% subplot(338),
% imagesc(surrSigmaC)
% set(gca,'TickDir','out','FontSize',14)
% title('Cov Cond. - Surr')
% axis square
% set(gcf,'color','w')
% 
% % Scatter plots comparing actual and surr data
% 
% 
% acT = reshape(targetSigmaT,1,[]);
% scT = reshape(surrSigmaT,1,[]);
% pT = polyfit(acT,scT,1);
% 
% subplot(333), hold on
% plot(acT,scT,'.k')
% set(gca,'TickDir','out','FontSize',14)
% axis square, box off
% xlabel('Cov Time - Actual'),ylabel('Cov Time - Surr')
% xl = xlim; yl = ylim;
% xlim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% ylim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% 
% xpT = [min(xl(1),yl(1)) max(xl(2),yl(2))];
% ypT = polyval(pT,xpT);
% plot(xpT,ypT,'k','linewidth',1)
% 
% 
% acN = reshape(targetSigmaN,1,[]);
% scN = reshape(surrSigmaN,1,[]);
% pN = polyfit(acN,scN,1);
% 
% subplot(336), hold on
% plot(acN,scN,'.k')
% set(gca,'TickDir','out','FontSize',14)
% axis square, box off
% xlabel('Cov Neurons - Actual'),ylabel('Cov Neurons - Surr')
% xl = xlim; yl = ylim;
% xlim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% ylim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% 
% xpN = [min(xl(1),yl(1)) max(xl(2),yl(2))];
% ypN = polyval(pN,xpN);
% plot(xpN,ypN,'k','linewidth',1)
% 
% 
% acC = reshape(targetSigmaC,1,[]);
% scC = reshape(surrSigmaC,1,[]);
% pC = polyfit(acC,scC,1);
% 
% subplot(339), hold on
% plot(acC,scC,'.k')
% set(gca,'TickDir','out','FontSize',14)
% xlabel('Cov Cond. - Actual'),ylabel('Cov Cond. - Surr')
% axis equal, axis square, box off
% xl = xlim; yl = ylim;
% xlim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% ylim([min(xl(1),yl(1)) max(xl(2),yl(2))])
% 
% xpC = [min(xl(1),yl(1)) max(xl(2),yl(2))];
% ypC = polyval(pC,xpC);
% plot(xpC,ypC,'k','linewidth',1)



% PLOT COVARIANCE MATRICES ACTUAL AND SURROGATE DATA, (same as previous but
% without scatter plot)

figure('units','normalized','outerposition',[0.2 0.1 .4 .8])
subplot(321),
imagesc(targetSigmaT)
set(gca,'TickDir','out','FontSize',14)
title('Cov Time - Actual')
axis square

subplot(323),
imagesc(targetSigmaN)
set(gca,'TickDir','out','FontSize',14)
title('Cov Neurons - Actual')
axis square

subplot(325),
imagesc(targetSigmaC)
set(gca,'TickDir','out','FontSize',14)
title('Cov Cond. - Actual')
axis square

subplot(322),
imagesc(surrSigmaT)
set(gca,'TickDir','out','FontSize',14)
ylabel('Cov Time')
axis square
title('Cov Time - Surr.')

subplot(324),
imagesc(surrSigmaN)
set(gca,'TickDir','out','FontSize',14)
title('Cov Neurons - Surr.')
axis square

subplot(326),
imagesc(surrSigmaC)
set(gca,'TickDir','out','FontSize',14)
title('Cov Cond. - Surr')
axis square
set(gcf,'color','w')



%  ========================================================================
%  ========================================================================
%%

% COMPARE TUNING

pars.tuning_conf_params = 'none'; %{'bootstrap',100,0.95};


% FOR THE ACTUAL DATA

% % get average firing rate in window -- FR x neuron
% fra = cell2mat(cellfun(@(x) sum(x,1)/(0.01*size(x,1)),{td.(pars.spiking_inputs{1})},'uni',0)');
% % get target directions - TARGET x 1
% dira = [td.target_direction];
% % fit tuning curves
% [tca,cba,ra] = regressTuningCurves(fra,dira',pars.tuning_conf_params);

% FOR THE TRIaL AVERAGED DATA
% get average firing rate in window -- FR x neuron
fra = cell2mat(cellfun(@(x) sum(x,1)/(0.01*size(x,1)),{tda.(pars.spiking_inputs{1})},'uni',0)');
% get target directions - TARGET x 1
dira = [tda.target_direction];
% fit tuning curves
[tca,cba,ra] = regressTuningCurves(fra,dira',pars.tuning_conf_params);


% FOR THE SURROGATE DATA -- I start with a tensor time x neuron x target

% get average firing rate in window -- FR x neuron
frsu = sum(surr_tensor,1)/(0.01*size(surr_tensor,1));
frsu = squeeze(frsu(1,:,:))';

% get target directions - TARGET x 1
dirsu = meta.targets;

% fit tuning curves
[tcsu,cbsu,rsu] = regressTuningCurves(frsu,dirsu',pars.tuning_conf_params);


% COMPARE TUNING CURVES
figure,
for i = 1:3
    subplot(1,3,i),hold on
    plot(tca(:,i),tcsu(:,i),'.k','markersize',16)
    xlabel('Actual data'),ylabel('Surrogate data')
    set(gcf,'color','w')
    set(gca,'TickDir','out','FontSize',14)
    if i == 1
        title('offset')
    elseif i == 2
        title('gain')
    else
        title('phase')
    end
end


% SEE IF TUNING CURVES LOOK SIMILIAr

% Normalize the firing rates by the max
nfra = zeros(size(fra,1),size(fra,2));
nfrsu = zeros(size(frsu,1),size(frsu,2));

for n = 1:size(fra,2)
    nfra(:,n) = fra(:,n)/max(fra(:,n));
end
for n = 1:size(frsu,2)
    nfrsu(:,n) = frsu(:,n)/max(frsu(:,n));
end


acol = parula(size(fra,2));


[~, ia] = sort(tca(:,3));
[~, isu] = sort(tcsu(:,3));

figure
subplot(121)
imagesc(nfra(:,ia)')
set(gca,'TickDir','out','FontSize',14)
xlabel('Target location')
ylabel('Neuron #')
set(gca,'XTick',[1 4 8],'XTickLabel',...
    {num2str(rad2deg(dira(1))) num2str(rad2deg(dira(4))) num2str(rad2deg(dira(8)))})
title('Actual data')
colorbar, caxis([0 1])
subplot(122)
imagesc(nfrsu(:,isu)')
set(gcf,'color','w')
set(gca,'TickDir','out','FontSize',14)
xlabel('Target location')
ylabel('Neuron #')
set(gca,'XTick',[1 4 8],'XTickLabel',...
    {num2str(rad2deg(dira(1))) num2str(rad2deg(dira(4))) num2str(rad2deg(dira(8)))})
title('Surrogate data')
colorbar, caxis([0 1])


% HISTOGRAM OF COSINE FITS
x_r = 0:0.05:1.05;
hra = histcounts(ra,x_r)/length(ra)*100;
hrsu = histcounts(rsu,x_r)/length(rsu)*100;

figure
hold on
ba = bar(x_r(1:end-1),hra,'histc');
bsu = bar(x_r(1:end-1),hrsu,'histc');
ba.FaceColor = 'r';
bsu.FaceColor = [.6 .6 .6];
alpha(ba,.5); alpha(bsu,.5)
set(gcf,'color','w')
set(gca,'TickDir','out','FontSize',14)
xlabel('R^2 cosine fits')
ylabel('Probability (%)')
xlim([0 1])
legend('Actual','Surr','Location','NorthWest'), legend boxoff
box off
title(['P = ' num2str(ranksum(ra,rsu)) ' (Wilcoxon)'])

% for n = 1:10:size(fra,2)
%     plot(dira,nfra(:,n),'color',acol(n,:));
% end
% set(gca,'TickDir','out','FontSize',14)
% subplot(122), hold on
% for n = 1:10:size(fra,2)
%     plot(dirsu,nfrsu(:,n),'color',acol(n,:))
% end
% set(gca,'TickDir','out','FontSize',14)
% set(gcf,'color','w')

% figure, hold on
% for n = 1:1:size(fra,2)
%     plot(nfra(:,n),nfrsu(:,n),'color',acol(n,:),'marker','.','markersize',16,'linestyle','none');
%     pause
% end
% axis square
% xlim([0 1]),ylim([0 1])