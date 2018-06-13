%
% QUICK & DIRTY CODE TO PLOT dPCA LATENT VARIABLES
%

% latent variable to plot
lvp = 5;

% data = dPCA_results struct
data = dr.lat_vars_mn;


% create time vector
bin_size = datasets{1}.stdata{1}.target{1}.bin_size;
time = 0:bin_size:bin_size*(size(data,4)-1);


% keep only the LVar to plot
data = squeeze(data(lvp,:,:,:));
dims = size(data);

% make colors for plotting
S    = size(data,1);
D    = size(data,2);
% create a number of different colormaps, for the conditions Si of each
% decision Di
colormap_list = {'summer','copper','autumn','winter','bone','pink','cool','hot'};
colors = zeros(S*D,3);
for i = 1:S
    eval(['tmp_color = ' colormap_list{i} '(D);'])
    for ii = 1:D
        colors((ii-1)*S+i,:) = tmp_color(ii,:);
    end
end

data = permute(data, [numel(dims) 1:numel(dims)-1]);
data = reshape(data, size(data,1), []);
data = data';

% this plots the 1st decision d1 for cond 1 (s1), for cond 2 (s2), ...
% cond n (sn), then 2nd decision d2 for cond 1 (s1), cond 2 (s2), ..
figure, hold on
for i = 1:S*D
    plot(time, data(i,:), 'LineWidth', 2, 'color', colors(i,:))
end
