%
% wrapper fcn to compare muscle synergy spaces across many pairs of tasks
% and sessions
%
%

function muscle_space_comp_all = compare_muscle_synergy_spaces_all( datasets, varargin )


% -------------------------------------------------------------------------
% read inputs
if nargin == 2
    % number of muscle synergies, necessary if doing PCA
    nbr_syn         = varargin{1};
end

% give error if muscle synergies were obtained with PCA, and haven't passed
% nbr_syn as input param
if strcmp(datasets{1}.dim_red_emg{1}.method,'pca')
    if ~exist('nbr_syn','var')
        error('PCA was used to get the muscle synergies, so you need to pass the nbr. of synergies as input param 2');
    end
elseif strcmp(datasets{1}.dim_red_emg{1}.method,'pca')
    if exist('nbr_syn','var')
        disp('the nbr. of synergies (param 2) will be disregarded for NMF')
    end
end


% -------------------------------------------------------------------------
% get nbr of datasets
nbr_dtsts           = length(datasets);


% do for all
for i = 1:nbr_dtsts
    % if used PCA
    if exist('nbr_syn','var')
        data{i}     = compare_muscle_synergy_spaces( datasets{i}, nbr_syn );
    % if used NMF
    else
        data{i}     = compare_muscle_synergy_spaces( datasets{i} );
    end
end


% return var
muscle_space_comp_all.data = data;



% -------------------------------------------------------------------------
% plot

% get nbr muscle synergies, if used NMF
if ~exist('nbr_syn','var')
    nbr_syn         = size(datasets{1}.dim_red_emg{1}.w,2);
end

colors              = distinguishable_colors(length(data));

figure, hold on
xlim([0 nbr_syn+1]),ylim([0 90])
set(gca,'Tickdir','out'),set(gca,'FontSize',16)
xlabel('dimension'),ylabel('angle (deg)')
for d = 1:length(data)
    plot( rad2deg(data{d}.princ_ang)', 'linewidth', 2, 'color', colors(d,:) );
end