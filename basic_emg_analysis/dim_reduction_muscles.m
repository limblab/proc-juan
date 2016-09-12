%
% Dimensionality reduction of EMGs using PCA or NNMF
%
% function dim_red_emg = dim_reduction_muscles( binned_data, varargin )
%
% Inputs (opt)          : [default]
%   binned_data         : binned_data struct or array of structs --can be
%                           cropped with 'cropped_bin_data.m' Alternatively
%                           this can be a matrix with EMGs (1 col is time)
%                           or a cell array of matrices of EMGs (again col
%                           1 is time, each task is a different cell)
%   (method)            : ['pca'] dim reduction method ('pca','nmf','none')
%   (chosen_emgs)       : ['all'] EMGs to be used for the analysis. 'all'
%                           will chose all
%   (labels)            : name of the task in each binned_data struct
%   (nbr_factors)       : nbr of factors we want to extract (only with NMF)
%   (plot_yn)           : [false] summary plot
%
% Outputs
%   dim_red_emg         : cell with fields:
%       w               : 'pca': eigenvectors in columns, 'nmf': weights in
%                           columns (muscles x nbr_factors)
%       eigen           : eigenvalues of w
%       scores          : result of applying w to smoothed_FR
%       t               : time axis for scores
%       chs             : EMG signals included in the analysis
%       method          : method used (stores input)
%
%


function dim_red_emg = dim_reduction_muscles( emg_data, varargin )


% ------------------------------------------------------------------------
% see what type of EMG input data has been passed, and reformat it to a
% struct called emg, with one databin and emgguide per input task

% see if it is a cell of either binned data structs or data matrices
if iscell( emg_data )

    emg(1:length(emg_data)) = struct('databin',[],'emgguide',[]);
    % see if it is a cell of binned data structs
    if isfield(emg_data{1},'emgguide')  
        for i = 1:length(emg_data)
            emg(i).databin = emg_data{i}.emgdatabin;
            emg(i).emgguide = emg_data{i}.emgguide;
            emg(i).timeframe = emg_data{i}.timeframe;
        end
    % or a cell of 2D matrices --note that column 1 is time
    else
        for i = 1:length(emg_data)
            emg(i).databin = emg_data{i}(:,2:end);
            emg(i).emgguide = cell(1,size(emg_data{i},2)-1);
            for ii = 1:(size(emg_data{i},2)-1)
                emg(i).emgguide{ii} = ['EMG' num2str(ii)];
            end
            emg(i).timeframe = emg_data{1}(:,1);
        end
    end
% or a single binned data struct
elseif isstruct( emg_data )
    emg.databin         = emg_data.emgdatabin;
    emg.emgguide        = emg_data.emgguide;
    emg.timeframe       = emg_data.timeframe;
% or a 2D array --note that column 1 is time
else
    % check if column 1 is time --if it is, get rid of it
    if length(unique(diff(emg_data(:,1)))) == 1
        emg.timeframe   = emg_data(:,1);
        emg_data(:,1)   = [];
    end
    emg.databin         = emg_data;
    for ii = 1:size(emg_data,2)
        emg.emgguide{ii} = ['EMG' num2str(ii)];
    end
end
clear emg_data;
        
        

% ------------------------------------------------------------------------
% get input params
if nargin == 1
    method              = 'pca';
    chosen_emgs         = 'all';
elseif nargin >= 2
    method              = varargin{1};
end
if nargin >= 3
    chosen_emgs         = varargin{2};
end
if nargin >= 4
    labels              = varargin{3};
end
if nargin >= 5
    nbr_factors         = varargin{4};
end
if nargin == 6
    plot_yn             = varargin{5};
end


% define var with the EMG channels that will be used
if isempty(chosen_emgs)
    chosen_emgs         = 'all';
end
if strcmp(chosen_emgs,'all')
    chosen_emgs         = 1:length(emg(1).emgguide);
end
if ~exist('chosen_emgs','var')
    chosen_emgs         = 1:length(emg(1).emgguide);
end
if ~exist('plot_yn','var')
    plot_yn             = false;
end


% clear labels if it's empty, so it won't break when plotting
if isempty(labels), clear labels; end


% nbr of BDFs
nbr_bdfs                = length(emg);
% nbr of EMGs --note that the input data may only comprise the chosen EMG
% channels, so check for that here
if size(emg(1).databin,2) == length(chosen_emgs);
    chosen_emgs         = 1:length(chosen_emgs);
end
nbr_emgs                = length(chosen_emgs);


% ------------------------------------------------------------------------
% Dim reduction
dim_red_emg             = cell(1,nbr_bdfs);
switch method
    % PCA
    case 'pca'
        for i = 1:nbr_bdfs
            [w_emg, scores_emg, eigen_emg] = pca(emg(i).databin(:,chosen_emgs));

            % store results
            dim_red_emg{i}.w        = w_emg;
            dim_red_emg{i}.eigen    = eigen_emg;
            dim_red_emg{i}.scores   = scores_emg;
            dim_red_emg{i}.t        = emg(i).timeframe;
            dim_red_emg{i}.chs      = chosen_emgs;
            dim_red_emg{i}.method   = 'pca';
            clear w_emg scores_emg eigen_emg
        end
    case 'nmf'
        if ~exist('nbr_factors','var')
            disp('you need to pass the number of factors for NMF');
        end
        for i = 1:nbr_bdfs
            [scores_emg, w_emg]     = nnmf(emg(i).databin(:,chosen_emgs),...
                                        nbr_factors);
                                    
            % store results
            dim_red_emg{i}.w        = w_emg';
            dim_red_emg{i}.scores   = scores_emg;
            if isfield(emg,'timeframe') % if not here will be filled by the calling fcn
                dim_red_emg{i}.t    = emg(i).timeframe;
            end
            dim_red_emg{i}.chs      = chosen_emgs;
            dim_red_emg{i}.method   = 'nmf';
            clear w_emg scores_emg
        end
    case 'none'
        for i = 1:nbr_bdfs
            
            % store results
            % -- the weights are just 1s because it is the raw EMGs
            dim_red_emg{i}.w        = ones(nbr_emgs);
            if size(emg(i).databin,2) ~= nbr_emgs
                dim_red_emg{i}.scores = emg(i).databin(:,chosen_emgs);
            else
                dim_red_emg{i}.scores = emg(i).databin;
            end
            dim_red_emg{i}.chs      = chosen_emgs;
            dim_red_emg{i}.method   = 'none';
        end
end



% ------------------------------------------------------------------------
% Return vars

% Turn a 1 element cell into an array
if length(dim_red_emg) == 1
    dim_red_emg                     = cell2mat(dim_red_emg);
end


% ------------------------------------------------------------------------
% Plot variance plots

if plot_yn
    
    switch method
        case 'pca'
            rows_plot               = floor(sqrt(nbr_bdfs));
            cols_plot               = ceil(sqrt(nbr_bdfs));
    
            f1h                     = figure;
            f2h                     = figure;

            % plots
            for i = 1:nbr_bdfs
                figure(f1h), subplot(rows_plot,cols_plot,i)
                bar(dim_red_emg{i}.eigen/sum(dim_red_emg{i}.eigen)),set(gca,'TickDir','out'),set(gca,'FontSize',14)
                xlim([0 nbr_emgs+1]),ylim([0 1])
                if exist('labels','var'), title(labels{i}); else title(['task #' num2str(i)]); end
                if rem(i-1,cols_plot) == 0
                    ylabel('norm. explained variance','FontSize',14)
                end
                if i >= ( (rows_plot-1)*cols_plot + 1 )
                    xlabel('component nbr.','FontSize',14)
                end

                figure(f2h), subplot(rows_plot,cols_plot,i)
                bar(cumsum(dim_red_emg{i}.eigen)/sum(dim_red_emg{i}.eigen))
                set(gca,'TickDir','out'),set(gca,'FontSize',14)
                xlim([0 nbr_emgs+1]),ylim([0 1])
                if exist('labels','var'), title(labels{i}); else title(['task #' num2str(i)]); end                
                if rem(i-1,cols_plot) == 0
                    ylabel('% norm. explained variance','FontSize',14)
                end
                if i >= ( (rows_plot-1)*cols_plot + 1 )
                    xlabel('component nbr.','FontSize',14)
                end
            end
            clear rows_plot cols_plot
            
        case 'nmf'
            colors_bars                 = parula(nbr_factors);
            
            figure('units','normalized','outerposition',[0 0 1 1])
            for i = 1:nbr_bdfs
                for ii = 1:nbr_factors
                    subplot(nbr_bdfs,nbr_factors,(i-1)*nbr_factors+ii)
                    bar(dim_red_emg{i}.w(ii,:),'FaceColor',colors_bars(ii,:))
                    set(gca,'TickDir','out'),set(gca,'FontSize',14)
                    xlim([0 nbr_emgs+1])
                    if ii == 1
                        ylabel(['weights --' labels{i}]);
                    end
                    if i == 1
                        title(['synergy ' num2str(ii)])
                    end
                    if i == nbr_bdfs
                        set(gca,'XTick',1:nbr_emgs)
                        set(gca,'XTickLabel',emg(1).emgguide(chosen_emgs))
                        set(gca,'XTickLabelRotation',45)
                    else
                        set(gca,'XTick',[]);
                    end
                end
            end
    end
end