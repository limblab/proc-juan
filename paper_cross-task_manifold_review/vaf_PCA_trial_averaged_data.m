%
% Do PCA on the trial-averaged data and compute the VAF
%


% what datasets to analyze
ds_use = 1:11;

% dimensionality of the manifold for which we want to compute VAF
mani_dim = 6;


for ds = 1:length(ds_use)
   
    nbr_units = length(datasets{ds_use(ds)}.neural_chs);
    nbr_tasks = length(datasets{ds_use(ds)}.labels);
    
    % pre-allocate matrices for storing the results -- define as NaNs
    % because due to trial-averaging some matrices seem not to be full rank
    % when doing PCA
    vaf_pca{ds}.norm_eigvals = nan(nbr_units,nbr_tasks);
	vaf_pca{ds}.cum_vaf = nan(nbr_units,nbr_tasks);
    
    % do for each task
    for t = 1:nbr_tasks
        
        % get trial-averaged data
        fr = datasets{ds_use(ds)}.stdata{t}.target{end}.neural_data.smoothed_fr_mn;
        
        % Do PCA
        [~, ~, eigvals] = pca(fr);
        
        % store VAF and cumsum(VAF)
        vaf_pca{ds}.monkey = datasets{ds_use(ds)}.monkey;
        vaf_pca{ds}.tasks = datasets{ds_use(ds)}.labels;
        vaf_pca{ds}.norm_eigvals(1:length(eigvals),t) = eigvals/sum(eigvals);
        vaf_pca{ds}.cum_vaf(1:length(eigvals),t) = cumsum(eigvals)/sum(eigvals);
        
        clear fr
    end
end


% -------------------------------------------------------------------------

% Get percentage VAF for a given manifold dimensionality

dims_explain_my_vaf = [];

for ds = 1:length(vaf_pca)
    
    for t = 1:length(vaf_pca{ds}.tasks)
        
        % and store
        dims_explain_my_vaf = [dims_explain_my_vaf, vaf_pca{ds}.cum_vaf(mani_dim,t)];
    end
end

% plot histogram
hist_bin = 0.05;
x_axis = 0:hist_bin:1+hist_bin;
[y_hist, ~] = histcounts(dims_explain_my_vaf,x_axis);

figure
bar(x_axis(1:end-1),y_hist/length(dims_explain_my_vaf)*100,'FaceColor',[.6 .6 .6])
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
ylabel('Percentage of tasks (%)')
xlabel(['Percentage neural variance explained with ' num2str(mani_dim) ' modes'])