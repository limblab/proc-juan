
clearvars -except datasets, close all;


% Manifold dimensionality
mani_dims_reach = 8;

% Params significance threshold
n_shuffles = 10000;
P_th = 0.001;


% what is what
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];


% preallocate matrices
all_PAs_reach = zeros(length(reach_ds),mani_dims_reach);
all_th_reach = zeros(length(reach_ds),mani_dims_reach);


%% ------------------------------------------------------------------------ 
% do for all reach-to-grasp tasks

for d = 1:length(reach_ds)
    
    % get eigenvenctors
    w1 = datasets{reach_ds(d)}.dim_red_FR{1}.w(:,1:mani_dims_reach);
    w2 = datasets{reach_ds(d)}.dim_red_FR{2}.w(:,1:mani_dims_reach);
    
    % compute PAs
    all_PAs_reach(d,:) = principal_angles(w1,w2);
    

    % compute significance threshold
    n_units = size(w1,1);
    [~, t_th] = empirical_principal_angle_distribution( n_units, mani_dims_reach, n_shuffles, P_th );
    
    all_th_reach(d,:) = t_th;
end




%% ------------------------------------------------------------------------ 

cols = parula(length(reach_ds)+1);

figure,hold on
for i = 1:length(reach_ds)
    plot(rad2deg(all_PAs_reach(i,:))','color',cols(i,:),'linewidth',2);
    plot(all_th_reach(i,:)','-.','color',cols(i,:),'linewidth',2);
end
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
xlabel('Neural mode'),ylabel('Principal angle (deg)');
xlim([0 mani_dims_reach]),ylim([0 90])

figure,hold on
for i = 1:length(reach_ds)
    plot((rad2deg(all_PAs_reach(i,:))./all_th_reach(i,:))','color',cols(i,:),'linewidth',2);
end
plot([0 mani_dims_reach],[1 1],'--','color',[.6 .6 .6],'linewidth',2)
set(gca,'TickDir','out','FontSize',14), box off,set(gcf,'color','w')
xlabel('Neural mode'),ylabel('Normalized Principal angle (deg)');
xlim([0 mani_dims_reach])