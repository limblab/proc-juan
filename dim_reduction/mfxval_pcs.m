%
% Do crossvalidation of PCs. The code works both cross-validating across
% time and neural channels.
%
%   stats = mfxval_pcs( bdf, neural_chs, xval_mode, xval_param, varargin )
%
% Inputs (opt)          : [default]
%   bdf                 : bdf struct 
%   neural_chs          : neural channels that define the complete data (in
%                           case the user wants to discard some)
%   xval_mode           : 'time' will do multifold cross-validation by
%                           dividing the data into n folds, then computing
%                           the correlation betwenn each of the folds and
%                           its corresponding portion of continuous data;
%                           'channels' will cross-validate by computing PC
%                           using m channels, projecting the data onto the
%                           time domain, and comparing these projections to 
%                           those obtained with the PCs computed using all
%                           the electrodes. This is repeated 'nbr_reps'
%                           times; 'subsets' will cross-validate btw
%                           subsets computed using m channels
%   xval_param          : if 'xval_mode' = 'time', length of the folds (s);
%                           if 'xval_mode' = 'channels', percentage of
%                           channels (0:1) that will be compared to all the
%                           channels; if 'xval_mode = subsets', percentage
%                           of channels (0:1) that will be used to compute
%                           the manifolds
%   (nbr_reps)          : [100] nbr_repetitions for xval_mode = 'channels', 
%                           
% Outputs:
%   stats               : canonical correlations and statistics comparing
%                           the PC projections with the complete set of
%                           electrodes vs. the subsampled sets
%
%
% Notes:
%   - There is some stuff in subsets that needs to be deleted
%

function stats = mfxval_pcs( bdf, neural_chs, xval_mode, xval_param, varargin )


% number of PCs that will be compared across the total and subsampled
% datasets --could be turn into a param
nbr_pcs_comp                = 20;
% bin size (s) --could be turn into a param
bin_size                    = 0.02;
% SD Gaussian kernel (s) --could be turn into a param
gauss_SD                    = 0.05;


% some preliminary definitions
% all_chs                     = double(arrayfun(@(x) x.id(1), bdf.units));
all_chs                     = 1:length(bdf.units);
discard_chs                 = setdiff(all_chs, neural_chs);

% -------------------------------------------------------------------------
% read input parameters
switch xval_mode
    case 'time'
        fold_length         = xval_param;
        if nargin == 5
            warning('input arg 4 will be ignored in time mode');
        end
    case {'channels','subsets'}
        perc_chs            = xval_param;
        if nargin == 5
            nbr_reps        = varargin{1};
        else
            nbr_reps        = 100;
        end
end

nbr_neural_chs              = numel(neural_chs);


% -------------------------------------------------------------------------
% do
switch xval_mode
    % --------------------------------------------------
    case 'time'
        
        % get nbr of folds
        if isfield(bdf,'emg')
            t_end           = bdf.emg.data(end,1);
            t_start         = bdf.emg.data(1,1);
        elseif isfield(bdf,'force')
            t_end           = bdf.force.data(end,1);
            t_start         = bdf.force.data(1,1);
        elseif isfield(bdf,'pos')
            t_end           = bdf.pos(end,1);
            t_start         = bdf.pos(1,1);
        else
            error('BDF does not have emg, force or pos fields');
        end
        nbr_folds           = floor((t_end-t_start)/fold_length);
        
        
        % preallocate matrices for scoring the stats
        R_CC                = zeros(nbr_pcs_comp,nbr_folds);
        R2_CC               = zeros(nbr_pcs_comp,nbr_folds);
        R2_raw              = zeros(nbr_pcs_comp,nbr_folds);
        
        % smooth firing rates and add them to the binned_data struct
        [smoothed_FR, binned_data] = gaussian_smoothing(bdf,'sqrt',bin_size,gauss_SD);
        dim_red_FR          = dim_reduction(smoothed_FR,'pca',discard_chs);
        
        
        disp(['Doing mfxval of the projections onto PCs over time...']);
        % do for each fold
        for f = 0:nbr_folds-1
            
            % 1. define the fold that will be left out
            data_out_start  = t_start + fold_length*f;
            data_out_end    = data_out_start + fold_length;
            
            % 2. get the training data window (will compute neural
            % primitives from here and test them on the synergies from the
            % entire data window
            [train_data, ~] = splitBinnedData( binned_data, ...
                                        data_out_start, data_out_end );
                                    
            % create a smoothed_FR array (with the time vector in col 1),
            % so it is compatible with dim_reduction
            train_smoothed_FR = [train_data.timeframe, train_data.smoothedspikerate];
                                    
            % 3. do PCA on train_data to compare its neural primitives to
            % those of the corresponding window obtained from 
            train_dim_red_FR = dim_reduction(train_smoothed_FR, ...
                                'pca',discard_chs);
                            
            % 4. get the neural primitives obtained from the entire dataset
            % that correspond to the primitives in train_dim_red_FR
            indx_scores_all = round(train_data.timeframe/bin_size);
            
            % ToDo: delete this dirty check
            if find( (dim_red_FR.t(indx_scores_all) - train_dim_red_FR.t)>0 )
                figure,plot(dim_red_FR.t(indx_scores_all) - train_dim_red_FR.t)
            end
            
            scores_all      = dim_red_FR.scores(indx_scores_all,1:nbr_pcs_comp);
            scores_train    = train_dim_red_FR.scores(:,1:nbr_pcs_comp);
            
            % 5. see how similiar the projections are 
            [~,~,R_CC(:,f+1),U,V,stats_CC{f+1}] = canoncorr(scores_all,scores_train);
            R2_CC(:,f+1)    = CalculateR2(U,V);
            
            R2_raw(:,f+1)   = CalculateR2(scores_all,scores_train);
            
            clear train_data train_smoothed_FR train_dim_red_FR scores_all scores_train;
        end
        
        % append statistics
        stats.R_CC          = R_CC;
        stats.mean_R_CC     = mean(R_CC,2);
        stats.std_R_CC      = std(R_CC,0,2);
        stats.stats_CC      = stats_CC;
        stats.R2_CC         = R2_CC;
        stats.mean_R2_CC    = mean(R2_CC,2);
        stats.std_R2_CC     = std(R2_CC,0,2);
        
        stats.R2_raw        = R2_raw;
        stats.mean_R2_raw   = mean(R2_raw,2);
        stats.std_R2_raw     = std(R2_raw,0,2);
  
        
    % --------------------------------------------------
    case 'channels'
        
        % preallocate matrix for storing random channels
        chs                 = zeros(nbr_reps,...
                                round(nbr_neural_chs*perc_chs));
        % preallocate matrix for storing CC correlations
        R_CC                = zeros(nbr_pcs_comp,nbr_reps);
        R2_CC               = zeros(nbr_pcs_comp,nbr_reps);
%         % preallocate matrix for storing corr with Procrustes 
%         R2_Proc             = zeros(nbr_pcs_comp,nbr_reps);
%         % stats_Proc          = zeros(nbr_pcs_comp,nbr_reps);

                        
        % smooth the FRs
        smoothed_FR         = gaussian_smoothing(bdf,'sqrt',bin_size,gauss_SD); %#ok<AGROW>                       
        % do PCA of all the selected channels
        dim_red_FR          = dim_reduction(smoothed_FR,'pca',discard_chs);
                        
        disp(['computing CCs of subsets of projections from ' num2str(perc_chs*100)...
            ' % of the electrodes...']);
        for i = 1:nbr_reps
            % randomly choose a subset of perc_chs
            chs(i,:)        = datasample(1:nbr_neural_chs,...
                            round(nbr_neural_chs*perc_chs),'Replace',false);
            discard_chs_this = setdiff(all_chs,sort(chs(i,:)));
            % and do pca
            dim_red_FR_xval = dim_reduction(smoothed_FR,'pca',discard_chs_this);
            % compare the dynamics of the projections using all the
            % channels, and the randonmly chosen subset of percentage of channels
            [~,~,R_CC(:,i),U,V,stats_CC{i}] = canoncorr(dim_red_FR.scores(:,1:nbr_pcs_comp),...
                            dim_red_FR_xval.scores(:,1:nbr_pcs_comp)); %#ok<AGROW>
            R2_CC(:,i)      = CalculateR2(U,V);
%             % run Procrustes algorithm 
%             [~,Z,transf]    = procrustes(dim_red_FR.scores(:,1:nbr_pcs_comp),...
%                                 dim_red_FR_xval.scores(:,1:nbr_pcs_comp));
%             R2_Proc(:,i)    = CalculateR2(dim_red_FR.scores(:,1:nbr_pcs_comp),Z);
%             % stats_Proc      = diag(P_mtrx_Proc);
            
        end
        
        % append statistics
        stats.R_CC          = R_CC;
        stats.mean_R_CC     = mean(R_CC,2);
        stats.std_R_CC      = std(R_CC,0,2);
        stats.stats_CC      = stats_CC;
        stats.R2_CC         = R2_CC;
        stats.mean_R2_CC    = mean(R2_CC,2);
        stats.std_R2_CC     = std(R2_CC,0,2);
        
%         stats.R_Proc        = R2_Proc;
%         % stats.stats_Proc    = stats_Proc;
%         stats.mean_R2_Proc  = mean(R2_Proc,2);
%         stats.std_R2_Proc   = std(R2_Proc,0,2);


    % --------------------------------------------------
    case 'subsets'

        % preallocate matrix for storing random channels
        chs_1               = zeros(nbr_reps,...
                                round(nbr_neural_chs*perc_chs));
        chs_2               = zeros(nbr_reps,...
                                round(nbr_neural_chs*perc_chs));
                            
        % preallocate matrix for storing CC correlations
        R_CC                = zeros(nbr_pcs_comp,nbr_reps);
        R2_CC               = zeros(nbr_pcs_comp,nbr_reps);

        % smooth the FRs
        smoothed_FR         = gaussian_smoothing2(bdf,'sqrt',bin_size,gauss_SD); %#ok<AGROW>                       
        
        disp(['Computing CCs btw pair of projections onto random subsets of ' ...
            num2str(perc_chs*100) ' % of the electrodes...']);
        
        
        for i = 1:nbr_reps
            % randomly choose the first subset of perc_chs
            chs_1(i,:)      = datasample(1:nbr_neural_chs,...
                            round(nbr_neural_chs*perc_chs),'Replace',false);
            chs_2(i,:)      = datasample(1:nbr_neural_chs,...
                            round(nbr_neural_chs*perc_chs),'Replace',false);
        
            % make sure they are not the same
            if chs_1 == chs_2, warning('two subsets of electrodes are identical!!!'); end
            
            % discard the channels when do not want (for making the random
            % subsets) ...
            discard_chs_1   = setdiff(all_chs,sort(chs_1(i,:)));
            discard_chs_2   = setdiff(all_chs,sort(chs_2(i,:)));
            % ...and do pca
            dim_red_FR_1    = dim_reduction(smoothed_FR,'pca',discard_chs_1);
            dim_red_FR_2    = dim_reduction(smoothed_FR,'pca',discard_chs_2);
            
            % ToDo: delete: shuffle dim_red_FR_2, it shouldn't change the
            % results if CCA is only based on the distributions
            rndm_indxs      = randperm(size(dim_red_FR_2.scores,1));
            for ii = 1:size(dim_red_FR_2.scores,2)
                dim_red_FR_2_shuffled.scores(:,ii) = dim_red_FR_2.scores(rndm_indxs,ii);
            end
            dim_red_FR_2_shuffled.scores = dim_red_FR_2.scores(rndm_indxs,:);
            
%            dim_red_FR_2.scores = dim_red_FR_2_shuffled.scores;
            % ToDo: delete up to here
            
            % compare the dynamics of the projections using all the
            % channels, and the randonmly chosen subset of percentage of channels
%             [~,~,R_CC(:,i),U,V,stats_CC{i}] = canoncorr(dim_red_FR_1.scores(:,1:nbr_pcs_comp),...
%                             dim_red_FR_2.scores(:,1:nbr_pcs_comp)); %#ok<AGROW>
            [A,B,R_CC(:,i),U,V,stats_CC{i}] = canoncorr(dim_red_FR_1.scores(:,1:nbr_pcs_comp),...
                            dim_red_FR_2.scores(:,1:nbr_pcs_comp)); %#ok<AGROW>
            R2_CC(:,i)      = CalculateR2(U,V);

           % append statistics
            stats.R_CC          = R_CC;
            stats.mean_R_CC     = mean(R_CC,2);
            stats.std_R_CC      = std(R_CC,0,2);
            stats.stats_CC      = stats_CC;
            stats.R2_CC         = R2_CC;
            stats.mean_R2_CC    = mean(R2_CC,2);
            stats.std_R2_CC     = std(R2_CC,0,2); 
            
            stats.A{i}          = A;
            stats.B{i}          = B;
            
            % ToDo: to delete
            [A_sh,B_sh,R_CC_sh(:,i),~,~,stats_CC_sh{i}] = canoncorr(dim_red_FR_1.scores(:,1:nbr_pcs_comp),...
                            dim_red_FR_2_shuffled.scores(:,1:nbr_pcs_comp)); %#ok<AGROW>
            stats.A_sh{i}       = A_sh;
            stats.B_sh{i}       = B_sh;
            stats.R_CC_sh{i}    = R_CC_sh;
            % ToDo: stop deleting here
        end
end

% -------------------------------------------------------------------------
% plot
switch xval_mode
    case 'time'
        figure,subplot(121),
        plot(R_CC,'color',[.6 .6 .6]), hold on
        plot(stats.mean_R_CC,'linewidth',4,'color','k')
        plot(stats.mean_R_CC+stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        plot(stats.mean_R_CC-stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        set(gca,'TickDir','out'),set(gca,'FontSize',14), ylim([0 1])
        xlabel('projection nbr.'),ylabel('canonical correlation')
        
        subplot(122),
        plot(R2_raw,'color',[1 .6 0]), hold on
        plot(stats.mean_R2_raw,'linewidth',4,'color','r')
        plot(stats.mean_R2_raw+stats.std_R2_raw,'linewidth',4,'color','r','linestyle','-.')
        plot(stats.mean_R2_raw-stats.std_R2_raw,'linewidth',4,'color','r','linestyle','-.')
        set(gca,'TickDir','out'),set(gca,'FontSize',14), ylim([0 1])
        xlabel('projection nbr.'),ylabel('R2 projections')
    
    case 'channels'
        % CCs of Canonical correlation
        figure,plot(R_CC,'color',[.6 .6 .6]), hold on
%         % R2 CCs
%         plot(R2_CC,'color',[1 .6 0]), hold on
        plot(stats.mean_R_CC,'linewidth',4,'color','k')
        plot(stats.mean_R_CC+stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        plot(stats.mean_R_CC-stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        
%         plot(stats.mean_R2_CC,'linewidth',4,'color','r')
%         plot(stats.mean_R2_CC+stats.std_R2_CC,'linewidth',4,'color','r','linestyle','-.')
%         plot(stats.mean_R2_CC-stats.std_R2_CC,'linewidth',4,'color','r','linestyle','-.')
        
        set(gca,'TickDir','out'),set(gca,'FontSize',14), ylim([0 1])
        xlabel('projection nbr.'),ylabel('canonical correlation')
        title(['CCs btw projections from ' num2str(perc_chs*100) ...
            ' % of the electrodes vs. all electrodes'])
    
    case 'subsets'
        figure,plot(R_CC,'color',[.6 .6 .6]), hold on
        plot(stats.mean_R_CC,'linewidth',4,'color','k')
        plot(stats.mean_R_CC+stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        plot(stats.mean_R_CC-stats.std_R_CC,'linewidth',4,'color','k','linestyle','-.')
        
        set(gca,'TickDir','out'),set(gca,'FontSize',14), ylim([0 1])
        xlabel('projection nbr.'),ylabel('canonical correlation')
        title(['CCs btw pair of projections onto random subsets of ' num2str(perc_chs*100)...
            ' % of the electrodes'])
end