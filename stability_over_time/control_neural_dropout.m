%
% Controls for stability over time paper
%


% Control 1. Neural population dynamics are robust to the specific channels
% used


% NEED TO HAVE RUN OCT_MANI_ALIGN2

reps = 100;
% percentages
percs = [.9 .8 .7 .6 .5];

t_s = 1;

[~,t_td] = getTDidx(master_td,{'date',sessions{t_s}});


% get spikes
sp = cat(1,t_td.M1_spikes);


% do

pool_cc = zeros(mani_dims(end),reps,length(percs));

for p = 1:length(percs)

    for r = 1:reps

        nu = size(sp,2);
        t_nu = floor(nu*percs(p));
        
        % Randomly split the neurons in two sets
        t_units = randperm(nu);
        units1 = t_units(1:t_nu);
        units2 = t_units((end-t_nu+1):end);
%         units1 = t_units(1:floor(nu/2));
%         units2 = t_units(floor(nu/2)+1:end);

        % do PCA
        [~, sc1] = pca(sp(:,units1));
        [~, sc2] = pca(sp(:,units2));

        % keep only as many PCs as the # of manifold dims
        sc1 = sc1(:,mani_dims);
        sc2 = sc2(:,mani_dims);

        [~,~,cc] = canoncorr(sc1,sc2);

        pool_cc(:,r,p) = cc;
    end
end



% Mean per percentage drop
mean_pool_cc = zeros(length(mani_dims),length(percs));
for p = 1:length(percs) 
   mean_pool_cc(:,p) =  mean(pool_cc(:,:,p),2);
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOOTSTRAPPING

cc_boots = zeros(mani_dims(end),reps);

for r = 1:reps
   
    [~, sc] = pca(sp);
    sc = sc(:,mani_dims);
    
    sc_shuf = sc(randperm(size(sc,1)),:);
    
    sc_shuf = smooth_data( sc_shuf, master_td(1).bin_size, 0.1);
    
    [~,~,cc_boots(:,r)] = canoncorr(sc,sc_shuf);
end

mean_boots_cc = mean(cc_boots,2);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS


% Legend
for p = 1:length(percs)
    lg{p} = [num2str(100-percs(p)*100) '% drop'];
end
lg{length(percs)+1} = 'shuffled in t';

% Raw data and means
cols_ctrl = parula(length(percs)+1);
figure, hold on
for p = 1:length(percs)
    for r = 1:reps
        if r == 1
            hp(p) = plot(pool_cc(:,r,p),'color',cols_ctrl(p,:));
        else
            plot(pool_cc(:,r,p),'color',cols_ctrl(p,:))
        end
    end
    plot(mean_pool_cc(:,p),'color',cols_ctrl(p,:),'linewidth',5);
end
for r = 1:reps
    if r == 1
        hp(length(percs)+1) = plot(cc_boots(:,r),'color',[.7 .7 .7]);
    else
        plot(cc_boots(:,r),'color',[.7 .7 .7])
    end
end
plot(mean_boots_cc,'color',[.7 .7 .7],'linewidth',5)
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+1])
xlabel('Neural mode')
ylabel('Similarity neural mode dynamics (CC)')
xlim([0 mani_dims(end)]); ylim([0 1])
hlg = [];
for p = 1:length(percs)
    hlg = [hlg, hp(p)];
end
legend(hlg,lg,'Location','West'), legend boxoff


% means only
figure, hold on
for p = 1:length(percs)
    plot(mean_pool_cc(:,p),'color',cols_ctrl(p,:),'linewidth',1.5);
end
plot(mean_boots_cc,'color',[.7 .7 .7],'linewidth',1.5)
set(gca,'TickDir','out','FontSize',14), box off,xlim([0 max(diff_days)+1])
xlabel('Neural mode')
ylabel('Similarity neural mode dynamics (CC)')
xlim([0 mani_dims(end)]); ylim([0 1])
legend(lg,'Location','West'), legend boxoff
