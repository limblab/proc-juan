

targs = unique([td.target_direction]);

plot_colors  =  distinguishable_colors(length(targs));

close all;
figure;

td_temp = td_all;
td2_temp = td_tune;
td3_temp = td_untune;
% 
% td_temp = trimTD(td_temp,'idx_movement_on',{'idx_movement_on',5});
% td2_temp = trimTD(td2_temp,'idx_movement_on',{'idx_movement_on',5});
% td3_temp = trimTD(td3_temp,'idx_movement_on',{'idx_movement_on',5});
td_temp = binTD(td_temp,'average');
td2_temp = binTD(td2_temp,'average');
td3_temp = binTD(td3_temp,'average');


ax(1) = subplot(1,3,1);
hold all;
for u = 1:length(targs)
    [~,temp] = getTDidx(td_temp,'target_direction',targs(u));
    temp = getSig(temp,{'M1_pca',1:3});
    plot3(temp(:,1),temp(:,2),temp(:,3),'.','MarkerSize',20,'Color',plot_colors(u,:));
end
axis square;

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
title('Full Population');

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');


ax(2) = subplot(1,3,2);
hold all;
for u = 1:length(targs)
    [~,temp] = getTDidx(td2_temp,'target_direction',targs(u));
    temp = getSig(temp,{'M1_pca',1:3});
    plot3(temp(:,1),temp(:,2),temp(:,3),'.','MarkerSize',20,'Color',plot_colors(u,:));
end
axis square;

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
title('Tuned Neurons');

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');


ax(3) = subplot(1,3,3);
hold all;
for u = 1:length(targs)
    [~,temp] = getTDidx(td3_temp,'target_direction',targs(u));
    temp = getSig(temp,{'M1_pca',1:3});
    plot3(temp(:,1),temp(:,2),temp(:,3),'.','MarkerSize',20,'Color',plot_colors(u,:));
end
axis square;

set(gca,'Box','off','TickDir','out','FontSize',14);
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
title('Untuned Neurons');

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');



linkprop(ax,'View');


if save_figs
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned_Clouds.fig']));
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned_Clouds.pdf']));
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned_Clouds.png']));
end



