%
% Generate empirical distributions of principal angles to repeat the
% analyses
%

% desired manifold dimensionalities
m = [8 10 15 20];

% parameters
P_orth = 0.001;
n_samples = 10000;

plot_yn = true;

% load paper datasets
load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');


% get space dimensionlities (number of units per dataset)
n = cellfun(@(x) numel(x.neural_chs), datasets);
clear datasets;

% do 
[dist_PAs, angle_non_orth] = empirical_principal_angle_distribution( n, m, n_samples, P_orth, plot_yn );


save(['/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/' ...
    'empirical principal angle distributions all datasets_revision']);