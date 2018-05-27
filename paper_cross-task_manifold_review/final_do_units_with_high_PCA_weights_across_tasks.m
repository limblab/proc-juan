%
% Code to try to address Reviewer 1 Question 5: "Are the units with high
% weights similar across tasks?"
%
%


% load the data
clearvars -except datasets;
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


mani_dims = 1:12;

% threshold high-weight units
th_high_w = 0.25; % in percentage


% initialize variables
all_high_w1 = [];
all_high_w2 = [];


this_ds = [1:3 7:9];% [4:6 10:11]; %[1:3 7:9];


for d = 1:length(this_ds)
    
    comb_tasks = nchoosek(1:length(datasets{this_ds(d)}.labels),2);
    % need to compare task a to b, and task b to a
    comb_tasks = [comb_tasks; fliplr(comb_tasks)];
    n_comb_tasks = size(comb_tasks,1);
    
    
    % do for all task combinations
    for c = 1:n_comb_tasks
        
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
        
        % retrieve weight distributions
        w1 = datasets{this_ds(d)}.dim_red_FR{t1}.w(:,mani_dims);
        w2 = datasets{this_ds(d)}.dim_red_FR{t2}.w(:,mani_dims);
        
        % -----------------------------------------------------------------
        % Find high weight threshold in task 1, and store high-weights in
        % task 1 and the weights of the same units onto the same mode in
        % task 2 
        
        % find threshold high weights in task 1
        all_w1 = reshape(w1,1,[]);
        high_w_th_t1 = prctile(abs(all_w1),(1-th_high_w)*100);
        
        % find position high weights, and retrieve weights
        % [r_h_w1,c_h_w1] = find(abs(w1)>high_w_th_t1);
        idx_h_w1 = find(abs(all_w1)>high_w_th_t1);
        h_w1 = w1(abs(w1)>high_w_th_t1);
        
        
        % find weights of those units in task 2
        all_w2 = reshape(w2,1,[]);
        h_w2 = all_w2(idx_h_w1)';
        
        % see percentage of high weight units in task 1 that have a high
        % weight in task 2
        high_w_th_t2 = prctile(abs(all_w2),(1-th_high_w)*100);
        idx_h_w2 = find(abs(all_w2)>high_w_th_t2);
        
        perc_h_w1_and_2 = sum(ismember(idx_h_w1,idx_h_w2))/length(idx_h_w1);
        
        
        % -----------------------------------------------------------------
        % Compare the weight of the units that have high-weights in task 1
        % with the weights that the same unit receives onto all modes from
        % task 2
        [r_h_w1, c_h_w1] = find(abs(w1)>high_w_th_t1);

        % retrieve the weights of each high-weight unit (for each mode) onto
        % all the modes from task 2
        h_w2_save = zeros(length(idx_h_w1),length(mani_dims));
        for i = 1:length(idx_h_w1)
            for j = 1:length(mani_dims)
                h_w2_save(i,j) = w2(r_h_w1(i),mani_dims(j));
            end
        end
        
        
        
        % -----------------------------------------------------------------
        % Prepare data for saving
        % repeat h_w1 as many times as modes, to compare this high weight
        % with the weight onto all modes
        h_w1_save = repmat(h_w1,1,12);
        h_w1_save = reshape(h_w1_save',1,[]);
        
        h_w2_save = reshape(h_w2_save',1,[]);
        
%         % plot for this comparison
%         figure,hold on,
%         plot(abs(h_w1),abs(h_w2),'.k')
%         yl = ylim;
%         plot([high_w_th_t1, high_w_th_t1],[0 yl(2)],'color',[.7 .7 .7],'linewidth',2)
%         plot([0 yl(2)],[high_w_th_t2, high_w_th_t2],'color',[.7 .7 .7],'linewidth',2)
%         xlim([0 yl(2)]),ylim([0 yl(2)])
%         set(gcf,'color','w')
%         set(gca,'TickDir','out','FontSize',14), box off
%         xlabel('Weight high-weight units on task 1')
%         ylabel('Weight same units on task 2')

        all_high_w1 = [all_high_w1, abs(h_w1_save)];
        all_high_w2 = [all_high_w2, abs(h_w2_save)];
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY COMPUTATIONS

all_high_w1 = double(all_high_w1);
all_high_w2 = double(all_high_w2);

xmax = ceil(max(all_high_w1)*10)/10;
ymax = ceil(max(all_high_w2)*10)/10;
max_ax = max(xmax,ymax);

[linfit, stats] = polyfit(all_high_w1,all_high_w2,1);
xfit = [0 max_ax];
yfit= polyval(linfit,xfit);


% subsample the points, because Corr will crash if they are too many
if length(all_high_w1) > 140000
    rnd_idx = round(rand([1 140000])*length(all_high_w1));
    rnd_idx = unique(rnd_idx);
else
    rnd_idx = 1:length(all_high_w1);
end
[r, P_r] = corr(all_high_w1(rnd_idx)',all_high_w2(rnd_idx)');


% Version 1
figure,hold on,
s = scatter(abs(all_high_w1),abs(all_high_w2),'filled','sizedata',10);
alpha(s,0.05);
s.CData = [0 0 0];
plot(xfit,yfit,'k','linewidth',1.5)
% yl = ylim;
% plot([high_w_th_t1, high_w_th_t1],[0 yl(2)],'color',[.7 .7 .7],'linewidth',2)
% plot([0 yl(2)],[high_w_th_t2, high_w_th_t2],'color',[.7 .7 .7],'linewidth',2)
% xlim([0 yl(2)]),ylim([0 yl(2)])
set(gcf,'color','w')
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Weight high-weight units on task 1')
ylabel('Weight same units on task 2')

% s = scatter(abs(all_high_w1),abs(all_high_w2),'filled','sizedata',10);
% alpha(s,0.05);
% s.CData = [1 0 0];
% plot(xfit,yfit,'r','linewidth',1.5)


% 
% % Version 2
% figure,hold on
% hist3([abs(all_high_w1),abs(all_high_w2)],'CdataMode','auto','FaceColor','interp','Nbins',[50 50])
% xlim([0.073 0.817]), ylim([0 0.74])
% plot3(xfit,yfit,[200 200],'--','color',[.6 .6 .6],'linewidth',3)
% colorbar
% view(2)
% set(gcf,'color','w')
% set(gca,'TickDir','out','FontSize',14), box off
% xlabel('Weight high-weight units on task 1')
% ylabel('Weight same units on task 2')
% 
