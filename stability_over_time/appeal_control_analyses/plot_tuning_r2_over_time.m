clear all;
close all;
clc;

root_dir = '/Users/mattperich/Dropbox/Research/Papers/M1 FES Grant documents/Juan 2018 - Stable latent activity/Results/Tuning stability/';

bins = 0:0.1:1;

n_rows = 1;
n_cols = 4;

f1 = figure('Position',[100 100 1000 400]);


%% Chewie
load(fullfile(root_dir,'NeuralTuningStability_Chewie_M1.mat'));

n_days = size(results.file_tuning,1);

hist_r = zeros(n_days, length(bins));
diff_days = zeros(1,n_days);

[all_r,all_r_labels,mean_r,std_r] = deal([]);
for iSess = 1:n_days
    r = mean(results.file_tuning{iSess,3},2);
    [n,x] = hist(r,bins);
    
    hist_r(iSess,:) = n;
    
    diff_days(iSess) = datenum(results.file_info(iSess).date,'mm-dd-yyyy') - datenum(results.file_info(1).date,'mm-dd-yyyy');
    
    all_r = [all_r; r];
    all_r_labels = [all_r_labels; diff_days(iSess)*ones(size(r))];
    mean_r = [mean_r; median(r)];
    std_r = [std_r; std(r)./sqrt(r)];
end

p  = anova1(all_r,all_r_labels,'off');
[~,~,~,~,s] = regress(all_r,[ones(size(all_r_labels)) all_r_labels]);


hist_r = 100*hist_r./repmat(sum(hist_r,2),1,size(hist_r,2));

figure(f1);
subplot(n_rows,n_cols,1);
imagesc(diff_days,bins,hist_r');
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YDir','normal');

xlabel('Days since first session');
ylabel('R^2');
title(['Chewie; r2 = ' num2str(s(1),3) '; p = ' num2str(s(3),3)]);




%% Chewie 2

load(fullfile(root_dir,'NeuralTuningStability_Chewie2_M1.mat'));

n_days = size(results.file_tuning,1);

[all_r,all_r_labels,mean_r,std_r] = deal([]);
hist_r = zeros(n_days, length(bins));
diff_days = zeros(1,n_days);
for iSess = 1:n_days
    r = mean(results.file_tuning{iSess,3},2);
    [n,x] = hist(r,bins);
    
    hist_r(iSess,:) = n;
    
    diff_days(iSess) = datenum(results.file_info(iSess).date,'mm-dd-yyyy') - datenum(results.file_info(1).date,'mm-dd-yyyy');
    all_r = [all_r; r];
    all_r_labels = [all_r_labels; diff_days(iSess)*ones(size(r))];
    mean_r = [mean_r; median(r)];
    std_r = [std_r; std(r)./sqrt(r)];
end

p  = anova1(all_r,all_r_labels,'off');
[~,~,~,~,s] = regress(all_r,[ones(size(all_r_labels)) all_r_labels]);

hist_r = 100*hist_r./repmat(sum(hist_r,2),1,size(hist_r,2));

figure(f1);
subplot(n_rows,n_cols,2);
imagesc(diff_days,bins,hist_r');
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YDir','normal');

xlabel('Days since first session');
ylabel('R^2');
title(['Chewie2; r2 = ' num2str(s(1),3) '; p = ' num2str(s(3),3)]);



%% Mihili

load(fullfile(root_dir,'NeuralTuningStability_Mihili_M1.mat'));

n_days = size(results.file_tuning,1);

[all_r,all_r_labels,mean_r,std_r] = deal([]);
hist_r = zeros(n_days, length(bins));
diff_days = zeros(1,n_days);
for iSess = 1:n_days
    r = mean(results.file_tuning{iSess,3},2);
    [n,x] = hist(r,bins);
    
    hist_r(iSess,:) = n;
    
    diff_days(iSess) = datenum(results.file_info(iSess).date,'mm-dd-yyyy') - datenum(results.file_info(1).date,'mm-dd-yyyy');
    all_r = [all_r; r];
    all_r_labels = [all_r_labels; diff_days(iSess)*ones(size(r))];
    mean_r = [mean_r; median(r)];
    std_r = [std_r; std(r)./sqrt(r)];
end

p  = anova1(all_r,all_r_labels,'off');
[~,~,~,~,s] = regress(all_r,[ones(size(all_r_labels)) all_r_labels]);

hist_r = 100*hist_r./repmat(sum(hist_r,2),1,size(hist_r,2));

figure(f1);
subplot(n_rows,n_cols,3);
imagesc(diff_days,bins,hist_r');
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YDir','normal');

xlabel('Days since first session');
ylabel('R^2');
title(['Mihili; r2 = ' num2str(s(1),3) '; p = ' num2str(s(3),3)]);



%% Jaco

load(fullfile(root_dir,'NeuralTuningStability_Jaco_M1.mat'));

n_days = size(results.file_tuning,1);

[all_r,all_r_labels,mean_r,std_r] = deal([]);
hist_r = zeros(n_days, length(bins));
diff_days = zeros(1,n_days);
for iSess = 1:n_days
    r = mean(results.file_tuning{iSess,3},2);
    [n,x] = hist(r,bins);
    
    hist_r(iSess,:) = n;
    
    diff_days(iSess) = datenum(results.file_info(iSess).date,'mm-dd-yyyy') - datenum(results.file_info(1).date,'mm-dd-yyyy');
    all_r = [all_r; r];
    all_r_labels = [all_r_labels; diff_days(iSess)*ones(size(r))];
    mean_r = [mean_r; median(r)];
    std_r = [std_r; std(r)./sqrt(r)];
end

p  = anova1(all_r,all_r_labels,'off');
[~,~,~,~,s] = regress(all_r,[ones(size(all_r_labels)) all_r_labels]);

hist_r = 100*hist_r./repmat(sum(hist_r,2),1,size(hist_r,2));

figure(f1);
subplot(n_rows,n_cols,4);
imagesc(diff_days,bins,hist_r');
set(gca,'Box','off','TickDir','out','FontSize',14);
axis tight;
set(gca,'YDir','normal');

xlabel('Days since first session');
ylabel('R^2');
title(['Jaco; r2 = ' num2str(s(1),3) '; p = ' num2str(s(3),3)]);



