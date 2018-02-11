%
% Compare EMG "manifolds" using principal angles
%


% Some definitions

mani_dim = 4; % based on the analysis of VAF

% wrist_sessions = [7 8 9 1:3]; % DARNED IT THE DATA FOR JACO ARE CRAP!!!
wrist_sessions = [7:9 1:3];
reach_sessions = [4:6 10:11];

sessions = [wrist_sessions, reach_sessions];

do_boots = true;
% bootstrapping params
boots_samples = 100000;
sign_boots = 0.001;

% -------------------------------------------------------------------------
% Do PCA of all the tasks --store data in pca_emg


% For the wrist tasks
for d = 1:length(sessions)
    
    n_tasks = length(datasets{sessions(d)}.labels);
    
    % Do PCA on the trial-related part of the data
    for t = 1:n_tasks

        % get concatenated single-trial EMGs
        emg = datasets{sessions(d)}.stdata{t}.target{end}.emg_data.conc_emg;    
        % do PCA
        w = pca(emg);
        % save in struct
        pca_emg{d}.w{t} = w;
    end

    % save in struct
    if ismember(d,wrist_sessions)
        pca_emg{d}.ds = 'wrist';
    else
        pca_emg{d}.ds = 'reach';
    end
end



% -------------------------------------------------------------------------
% Principal angles between the EMG manifolds from all tasks in each session

for d = 1:length(pca_emg)
   
    % get all task pairs
    combs = nchoosek(1:length(pca_emg{d}.w),2);
    
    % preallocate matrix
    PA{d} = zeros(size(combs,1),mani_dim);
    
    % compute principal angles between their manifolds
    for c = 1:size(combs,1)
        
        PA{d}(c,:) = principal_angles( pca_emg{d}.w{combs(c,1)}(:,1:mani_dim), ...
                        pca_emg{d}.w{combs(c,2)}(:,1:mani_dim) );
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do bootstrapping to establish significance

if do_boots
    
    % start parallel computing
    gcp;
    
    % get space dimensionality across all datasets
    emg_chs = cellfun(@(x) length(x.chosen_emgs), datasets);
    space_dim = sort(unique(emg_chs(sessions)));
    
    % do!
    [dist_angles, angle_non_orth] = empirical_principal_angle_distribution( space_dim, mani_dim, boots_samples, sign_boots, true );
    
    % squeeze the matrix because we are always using the same manifold
    % dimensionality
    angle_non_orth = squeeze(angle_non_orth(:,1,:));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the largest PA that is significantly small

% pre-allocate matrix to save the results
largest_align_dim = [];

for d = 1:length(sessions)
   
    t_th = angle_non_orth(:,space_dim==length(datasets{sessions(d)}.chosen_emgs))';

    % do not consider 
    
    % find largest non-orthogonal mode and store it
    for c = 1:size(PA{d},1)
        
        non_orth_mode = find(rad2deg(PA{d}(c,:)) > t_th, 1) -1;
        % if it's empty it means all EMG manifold dimensions are
        % well-aligned
        if isempty(non_orth_mode)
            non_orth_mode = mani_dim;
        end
        
        largest_align_dim = [largest_align_dim, non_orth_mode];
    end
end

% do a histogram
x_axis = 0:space_dim+1;
y_hist = histcounts(largest_align_dim,x_axis);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS

% Replicate figure 1 in the paper
ds_wrist = find(sessions==wrist_sessions(2));
ds_reach = find(sessions==reach_sessions(end));

figure,
subplot(131), hold on
plot(rad2deg(PA{ds_reach})','linewidth',1.5)
plot(angle_non_orth(:,space_dim==length(datasets{ds_reach}.chosen_emgs)),...
    'linewidth',3,'linestyle','-.','color',[.5 .5 .5])
ylim([0 90]),xlim([0 mani_dim])
ylabel('Principal angle (deg)')
xlabel('EMG mode'), title('Reach-to-grasp')
set(gca,'TickDir','out','FontSize',14), box off
subplot(132), hold on
plot(rad2deg(PA{ds_wrist})','linewidth',1.5)
plot(angle_non_orth(:,space_dim==length(datasets{ds_wrist}.chosen_emgs)),...
    'linewidth',3,'linestyle','-.','color',[.5 .5 .5])
ylim([0 90]),xlim([0 mani_dim])
xlabel('EMG mode'), title('Wrist')
set(gca,'TickDir','out','FontSize',14), box off
subplot(133)
bar(x_axis(1:end-1),y_hist,'FaceColor','k')
set(gca,'TickDir','out','FontSize',14), box off
yl = ylim; xl = xlim; text(xl(2)-2,yl(2)-1,['n=' num2str(length(largest_align_dim))],'FontSize',14)
xlabel(['Highest non-orth EMG mode (P<' num2str(sign_boots)]),ylabel('Counts')
