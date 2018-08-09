

figure('Position',[100 100 800 400]);

subplot(1,3,1);
hold all;
temp = confusionmat(clas_results.dir_within(1,:),clas_results.pred_within(1,:));
temp = temp./repmat( sum(temp,2),1,size(temp,2) );
imagesc(temp);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
set(gca,'XTick',1:8,'YTick',1:8);
% set(gca,'CLim',[0 23]);
xlabel('True target');
ylabel('Predicted target');
axis tight;
axis square;
colormap(brewermap(23,'*Blues'));
title('Within');
colorbar

subplot(1,3,2);
hold all;
temp = confusionmat(clas_results.dir_across(1,:),clas_results.pred_across(1,:));
temp = temp./repmat( sum(temp,2),1,size(temp,2) );
imagesc(temp);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
set(gca,'XTick',1:8,'YTick',1:8);
xlabel('True target');
ylabel('Predicted target');
axis tight;
axis square;
title('Across');
colorbar

subplot(1,3,3);
hold all;
temp = confusionmat(clas_spike_results.dir_across_spike(1,:),clas_spike_results.pred_across_spike(1,:));
temp = temp./repmat( sum(temp,2),1,size(temp,2) );
imagesc(temp);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
set(gca,'XTick',1:8,'YTick',1:8);
xlabel('True target');
ylabel('Predicted target');
axis tight;
axis square;
title('Spike')
colorbar

