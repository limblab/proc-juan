%
% Compare cross-task and within-task principal angles
%


% Load the data
load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
load('/Users/juangallego/Documents/Publications/2017 - Multi task manifold + some stuff not on DB/Raw data for figs - 2017-06-20/Results_manifold_notebook_v1_20170620.mat');

% dataset and tasks to compare
d = 7;
t1 = 1;
t2 = 2;

% Parameters
reps = 1000; % repetitions for shuffling
mani_dim = 12;

% -------------------------------------------------------------------------
% Randomly split the trials in tasks 1 and 2 in two 

sfr1 = datasets{d}.stdata{t1}.target{end}.neural_data.smoothed_fr;
sfr2 = datasets{d}.stdata{t2}.target{end}.neural_data.smoothed_fr;

n_trials1 = size(sfr1,3);
n_trials2 = size(sfr2,3);


% Do
PA1 = zeros(mani_dim,reps);
PA2 = zeros(mani_dim,reps);


for s = 1:reps    
    
    % take two random subsets of trials
    r_trials1 = randperm(n_trials1,floor(n_trials1/2));
    o_trials1 = setdiff(1:n_trials1,r_trials1);
    o_trials1 = o_trials1(1:length(r_trials1));
    
    % Prepare matrix to have two sets of concatenated trials
    r_sfr1 = sfr1(:,:,r_trials1);
    o_sfr1 = sfr1(:,:,o_trials1);
    
    pr_sfr1 = permute(r_sfr1,[1 3 2]);
    po_sfr1 = permute(o_sfr1,[1 3 2]);
    
    cr_sfr1 = reshape(pr_sfr1,[],size(r_sfr1,2));
    co_sfr1 = reshape(po_sfr1,[],size(o_sfr1,2));
    
    % Do PCA
    w_r = pca(cr_sfr1);
    w_o = pca(co_sfr1);
    
    % compute principal angles between manifolds
    tPA = principal_angles(w_r(:,1:mani_dim),w_o(:,1:mani_dim));
    
    % and store them
    PA1(:,s) = tPA;
end


% Load or Compute cross-task PA
if exist('angle_results','var')
    
    combs = nchoosek(1:length(datasets{d}.labels),2);
    [~, b] = ismember(combs,[t1 t2],'rows');
    cPA = angle_results.data{d}.princ_angles.angles(find(b),:);
    disp('using saved cross-task PA')
else
    cPA = principal_angles(datasets{d}.dim_red_FR{t1}.w(:,1:mani_dim),...
                        datasets{d}.dim_red_FR{t2}.w(:,1:mani_dim));
end

% plot and compare to cross-task PA
figure,hold on
plot(rad2deg(cPA),'k','linewidth',2)
plot(rad2deg(mean(PA1,2)),'r','linewidth',2)
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Neural mode'),ylabel('Principal Angle (deg')
if exist('angle_results','var')
    plot(rad2deg(angle_results.data{d}.princ_angles.angle_orth),'color',[.6 .6 .6],...
        'linewidth',2)
    legend('cross-task PA','within-task control','random control','Location','NorthWest')
else
    legend('cross-task PA','within-task control','Location','NorthWest')
end
legend boxoff
title([datasets{d}.monkey ' ' datasets{d}.date(1:end-4) ' - ' datasets{d}.labels{1} ' vs ' datasets{d}.labels{2}])