%
% Trying to understand the dependence of the principal angles on the
% analysis window
%
% -- This analysis is based on my observation that taking the entire
% reach-to-grasp trials seems to make the PAs larger than if we take the
% several first hundred ms, and on Surya's idea that Neural Task Complexity
% increases linearly with the analysis window, which may in turn change
% dimensionality (?)
%


% Parameters
mani_dims = 12;
nbr_bins = [20 25 30 35]; % take the first N bins
ds = [1:3 7:9]; % What datasets; WRIST: ds = [1:3 7:9]; REACH-TO-GRASP: [4:6 10:11]

% load for the cross-task manifold paper
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do for all sessions and pairs of tasks
for d = 1:length(ds)
    
    % get combinations of tasks for this sessions
    combs_t = nchoosek(1:length(datasets{ds(d)}.labels),2);
    n_combs_t = size(combs_t,1);

    
    % preallocate cell for 
    for c = 1:n_combs_t
        
        all_PAs{d}.task_pair{c}.data = zeros(mani_dims,length(nbr_bins));
    end
        
    
    % do for all task combinatinos
    for c = 1:n_combs_t

        
        rsfr1 = datasets{ds(d)}.stdata{combs_t(c,1)}.target{end}.neural_data.smoothed_fr;
        rsfr2 = datasets{ds(d)}.stdata{combs_t(c,2)}.target{end}.neural_data.smoothed_fr;

        % -----------------------------------------------------------------
        % do for all bin sizes
        for b = 1:length(nbr_bins)
            
            % take the desired analysis windows 
            sfr1 = rsfr1(1:nbr_bins(b),:,:);
            sfr2 = rsfr2(1:nbr_bins(b),:,:);

            % Prepare the data to do PCA: create a (time x trials) x neurons
            % matrix
            sfr1 = permute(sfr1,[1 3 2]);
            sfr2 = permute(sfr2,[1 3 2]);

            sfr1 = reshape(sfr1,[],size(sfr1,3));
            sfr2 = reshape(sfr2,[],size(sfr2,3));
            
            % -------------------------------------------------------------
            % Do PCA and compute the Principal Angles
            
            w1 = pca(sfr1);
            w2 = pca(sfr2);
            
            w1 = w1(:,1:mani_dims);
            w2 = w2(:,1:mani_dims);
            
            PA = principal_angles(w1,w2);
            
            % Save the data
            all_PAs{d}.task_pair{c}.data(:,b) = PA;
        end
        
        clear rsfr* sfr*
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary analyses
