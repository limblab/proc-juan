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
% compute pairwise angles

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


% Compute distributions of absolute value of dot products (so angles lie
% within 0-90)
h1_bin = 0.05;
x_h1 = 0:h1_bin:1+h1_bin;
y_h1 = histcounts(abs(reshape(dot_prods,1,[])),x_h1)/numel(dot_prods)*100;

figure,
b1 = bar(x_h1(1:end-1),y_h1,'histc');
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 1])
set(b1,'FaceColor',[.8 .8 .8])
ylabel('Pairwise dPC comparisons (%)')
xlabel('|Dot product|')


% For the angles themselves
h2_bin = 5;
x_h2 = 0:h2_bin:90+h2_bin;
cos_dot_prods = rad2deg( acos( abs(reshape(dot_prods,1,[])) ) );
y_h2 = histcounts(cos_dot_prods,x_h2)/numel(cos_dot_prods)*100;

mn_cos = mean(cos_dot_prods);
sd_cos = std(cos_dot_prods);

yl2 = ceil(max(y_h2)*1.2);
ystats = (max(y_h2)+yl2)/2;

c2 = [.8 .8 .8];

figure,hold on
b2 = bar(x_h2(1:end-1),y_h2,'histc');
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 90])
set(b2,'FaceColor',c2)
ylabel('Pairwise dPC comparisons (%)')
xlabel('Angle between dPCs')
ylim([0 yl2])
plot(mn_cos,ystats,'o','markersize',12,'color',c2)
plot([mn_cos-sd_cos, mn_cos+sd_cos],[ystats ystats],'linewidth',1.5,'color',c2);
text(10,yl2*.85,['n=' num2str(numel(cos_dot_prods))],'Fontsize',14)


% For the example dataset
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