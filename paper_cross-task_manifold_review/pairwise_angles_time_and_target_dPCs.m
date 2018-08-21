%
% Summarize pairwise angles between dPCs
%

dPCA_ds = [1:3 7:9];
mani_dims = 12;

% load all the data
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% Do dPCA if not in the ws
if ~exist('dPCA_results','var')
    for d = 1:length(dPCA_ds)
        dPCA_results{d} = call_dPCA( datasets{dPCA_ds(d)}.stdata, mani_dims, false );
    end
end


% -------------------------------------------------------------------------
% Some definitions

time_marg = find(strcmpi(dPCA_results{1}.marg_names,'time'));
target_marg = find(strcmpi(dPCA_results{1}.marg_names,'target'));

% Alignment threshold (P<0.001 being less differet than random)
% dot prod > 3.3*N^(-.5)
align_th = zeros(1,length(dPCA_ds));
for d = 1:length(dPCA_ds)
    % space dimensionality
    N = size(dPCA_results{d}.W,1);
    % smallest dot product by chance; and smallest dot product by chance
    dot_th = 3.3/sqrt(N);
    align_th(d) = acosd(dot_th);
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% compute pairwise angles between time-related and target-related 

ctr = 1;

for d = 1:length(dPCA_ds)
   
    % find time-related modes
    time_modes = find(dPCA_results{d}.which_marg == time_marg);
    target_modes = find(dPCA_results{d}.which_marg == target_marg);
    
    % compute angle between each time-related mode and each target-related mode
    for tm = 1:length(time_modes)
        
        for tgtm = 1:length(target_modes)
            
            t_time_mode = time_modes(tm);
            t_target_mode = target_modes(tgtm);
            
            utime = dPCA_results{d}.W(:,t_time_mode) / norm(dPCA_results{d}.W(:,t_time_mode));
            utarget = dPCA_results{d}.W(:,t_target_mode) / norm(dPCA_results{d}.W(:,t_target_mode));
            
            dot_prods_time_target{d}(ctr) = dot(utime,utarget);
            
            ctr = ctr + 1;
        end
    end
    
    ctr = 1;
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% compute pairwise angles between all dPCs

ctr = 1;

for d = 1:length(dPCA_ds)
    
    % compute angle between dim m and dim m+1, m+2, ...
    for m = 1:mani_dims
        
        for m2 = m+1:mani_dims
            
            u1 = dPCA_results{d}.W(:,m) / norm(dPCA_results{d}.W(:,m));
            u2 = dPCA_results{d}.W(:,m2) / norm(dPCA_results{d}.W(:,m2));
            
            dot_prods(d,ctr) = dot(u1,u2);
            ctr = ctr + 1;
        end
    end
    
    ctr = 1;
end


% -------------------------------------------------------------------------
% plot


% % Compute distributions of absolute value of dot products (so angles lie
% % within 0-90)
% h1_bin = 0.05;
% x_h1 = 0:h1_bin:1+h1_bin;
% y_h1 = histcounts(abs(reshape(dot_prods,1,[])),x_h1)/numel(dot_prods)*100;
% 
% figure,
% b1 = bar(x_h1(1:end-1),y_h1,'histc');
% set(gcf, 'color', [1 1 1])
% set(gca,'FontSize',14,'TickDir','out'); box off
% xlim([0 1])
% set(b1,'FaceColor',[.8 .8 .8])
% ylabel('Pairwise dPC comparisons (%)')
% xlabel('|Dot product|')


% For the angles themselves
h2_bin = 5;
x_h2 = 0:h2_bin:90+h2_bin;
acos_dot_prods = rad2deg( acos( abs(reshape(dot_prods,1,[])) ) );
y_h2 = histcounts(acos_dot_prods,x_h2)/numel(acos_dot_prods)*100;

acos_dot_timetarget = rad2deg(acos(abs(cell2mat(dot_prods_time_target))));
y_htt = histcounts(acos_dot_timetarget,x_h2)/numel(acos_dot_timetarget)*100;


mn_cos = mean(acos_dot_prods);
sd_cos = std(acos_dot_prods);

mn_cos_timtarget = mean(acos_dot_timetarget);
sd_cos_timtarget = std(acos_dot_timetarget);


yl2 = ceil(max(y_h2)*1.2);
ystats = (max(y_h2)+yl2)/2;

c2 = [.8 .8 .8];


% % ANGLES BETWEEN ALL PAIRS OF dPCS
% figure,hold on
% b2 = bar(x_h2(1:end-1),y_h2,'histc');
% set(b2,'FaceColor',c2)
% set(gcf, 'color', [1 1 1])
% set(gca,'FontSize',14,'TickDir','out'); box off
% xlim([0 90])
% ylabel('Pairwise dPC comparisons (%)')
% xlabel('Angle between dPCs')
% ylim([0 yl2])
% plot(mn_cos,ystats,'o','markersize',12,'color',c2)
% plot([mn_cos-sd_cos, mn_cos+sd_cos],[ystats ystats],'linewidth',1.5,'color',c2);
% text(10,yl2*.85,['n=' num2str(numel(acos_dot_prods))],'Fontsize',14)


% ANGLE BETWEEN TIME- AND TARGET-RELATED DPCS
figure,hold on
b3 = bar(x_h2(1:end-1),y_htt,'histc');
set(b3,'FaceColor','b')
alpha(b3,.5);
plot([max(align_th) max(align_th)],[0 yl2],'r','linewidth',2)
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 90])
ylabel('Pairwise comparisons time- vs target-dPCs (%)')
xlabel('Angle between dPCs')
ylim([0 yl2])
plot(mn_cos_timtarget,ystats+1,'o','markersize',12,'color','b')
plot([mn_cos_timtarget-sd_cos_timtarget, mn_cos_timtarget+sd_cos_timtarget],[ystats+1 ystats+1],'linewidth',1.5,'color','b');
text(10,yl2*.85,['n=' num2str(numel(acos_dot_timetarget))],'Fontsize',14,'color','b')
legend('pairwise angles','largest non-orth th','Location','West'),legend boxoff

% count the number of pairwise angles btw time- and target-related dPCs
% that are below the non-orth threshold 
n_non_orth = zeros(1,length(dPCA_results));
n_angles = zeros(1,length(dPCA_results));
for d = 1:length(dPCA_results)
   
    n_non_orth(d) = sum(rad2deg(acos(abs(dot_prods_time_target{d}))) < align_th(d));
    n_angles(d) = length(dot_prods_time_target{d});
end
perc_non_orth = sum(n_non_orth)/sum(n_angles)*100;
disp(['% non-orthogonal time-related and target-related dPCs: ' num2str(perc_non_orth,2)]);


% COMPARISON ANGLES BETWEEN ALL PAIRS OF VECTORS AND ANGLES BETWEEN TIME
% AND TARGET -RELATED MODES
% -- Note that All also includes time vs target
figure,hold on
b4 = bar(x_h2(1:end-1),y_h2,'histc');
set(b4,'FaceColor',c2)
b5 = bar(x_h2(1:end-1),y_htt,'histc');
set(b5,'FaceColor','b')
alpha(b5,.5);
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 90])
ylabel('Pairwise dPC comparisons (%)')
xlabel('Angle between dPCs')
ylim([0 yl2])
plot(mn_cos,ystats,'o','markersize',12,'color',c2)
plot([mn_cos-sd_cos, mn_cos+sd_cos],[ystats ystats],'linewidth',1.5,'color',c2);
plot(mn_cos_timtarget,ystats+1,'o','markersize',12,'color','b')
plot([mn_cos_timtarget-sd_cos_timtarget, mn_cos_timtarget+sd_cos_timtarget],[ystats+1 ystats+1],'linewidth',1.5,'color','b');
text(10,yl2*.85,['n=' num2str(numel(acos_dot_prods))],'Fontsize',14)
text(10,yl2*.75,['n=' num2str(numel(acos_dot_timetarget))],'Fontsize',14,'color','b')
legend('all dPC pairs','time vs target','Location','West'), legend boxoff


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% For an example dataset
ex_ds = 5;

% Draw random distribution of pairwise angles -the distirbution of dot
% products is a normal
n_ex = size(dPCA_results{ex_ds}.W,1);
x_rand_dot = [-3.3:.01:3.3];
rand_dot_ex = normpdf(x_rand_dot,0,1/sqrt(n_ex));
rand_angles_ex = abs(acosd(abs(rand_dot_ex))); 

% get distribution PAs for example dataset
h3_bin = h2_bin;
x_h3 = 0:h3_bin:90+h3_bin;
cos_dot_prods_ex = rad2deg(acos(abs(reshape(dot_prods(ex_ds,:),1,[]))));
y_h3 = histcounts(cos_dot_prods_ex,x_h3)/numel(cos_dot_prods_ex)*100;

nonorth = rad2deg(3.3/sqrt(size(dPCA_results{ex_ds}.W,1)));

mn_cos_ex = mean(cos_dot_prods_ex);
sd_cos_ex = std(cos_dot_prods_ex);

% get distribution random angles
x_rnd = 0:1:91;
y_rnd = histcounts(rand_angles_ex,x_rnd)/numel(rand_angles_ex)*100;

yl3 = ceil(max(y_h3)*1.2);
ystats3 = (yl3+max(y_h3))/2;

figure, hold on
b3 = bar(x_h3(1:end-1),y_h3,'histc');
plot(x_rnd(1:end-1),y_rnd,'r','linewidth',1.5)
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 90])
set(b3,'FaceColor',c2)
ylabel('Pairwise dPC comparisons (%)')
xlabel('Angle between dPCs')
ylim([0 yl3])
plot(mn_cos_ex,ystats3,'o','markersize',12,'color',c2)
plot([mn_cos_ex-sd_cos_ex, mn_cos_ex+sd_cos_ex],[ystats3 ystats3],'linewidth',1.5,'color',c2);
text(10,yl3*.85,['n=' num2str(numel(cos_dot_prods_ex))],'Fontsize',14)
title([datasets{dPCA_ds(ex_ds)}.monkey ' - ' datasets{dPCA_ds(ex_ds)}.date])
legend('angles','random','Location','West'), legend boxoff