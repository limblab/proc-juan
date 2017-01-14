%
% Compare EMG fits using decoders that take different sets
% (marginalizations) of dPCs as inputs
%

clearvars -except dPCA_data datasets

% check that the parallel pool is running, otherwise start it
gcp;


% ---------------------------------------------------------------
% Some definitions

% crossvalidate the predicitons?
xval_yn     = true;
% manifold dimensionality
mani_dim    = 12;
% datasets that will be analyzed
D           = [1:3 7:9]; % Jaco and Jango 1D / 2D datasets


% marginalization legends for plotting fits
marg_leg    = {'task','target','time','task/tgt int'};

% combinations of marginalizations that will be compared
marg_set    = {1, 2, 3, 4, [2 4], [1 3], [2 3], [1 2 4], [1 2 3], [1 3 4], [1 2 3 4]};



% ---------------------------------------------------------------
% 0. Do dPCA of the neural data, if the data are not avaialble

if ~exist('dPCA_data','var')
    for d = 1:length(D)
        dPCA_data{d} = call_dPCA( datasets{D(d)}.stdata, mani_dim );
    end
end



% ---------------------------------------------------------------
% 1. Build the decoders 

% Build decoders taking the leading 'dpcs_marg' dPCs per marginalization as
% predictors

% number of dPCs per marginalization to use
dpcs_marg   = 2;

for d = 1:length(D)
    disp(['building decoders with chosen marginalizations for dataset ', num2str(D(d))]);
    for m = 1:length(marg_set)
        dPCA_fit{d,m} = build_dPC_decoder( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'margs',            marg_set{m}, ...
                            'nbr_dpcs_p_marg',  dpcs_marg, ... 
                            'xval_yn',          xval_yn, ...
                            'fold_size',        60, ...
                            'smooth_spikes',    true, ...
                            'dec_offset',       true, ...
                            'poly_order',       2, ...
                            'return_preds',     false, ...
                            'plot_yn',          false );
    end
end
disp('...');


% --------------------
% Build decoders taking random combinations of dPCS as predictors

% number of random groups
nbr_rnd_iter = 10;

% get minimum and maximum number of input marginalizations
max_nbr_marg = max(cellfun(@(x) numel(x), marg_set));
min_nbr_marg = min(cellfun(@(x) numel(x), marg_set));

% define random combinations of dPCs to use as predictors --will be the
% same for all the datasets
for m = min_nbr_marg:max_nbr_marg
    for i = 1:nbr_rnd_iter
        rand_dpcs{m,i} = datasample(1:mani_dim,m*dpcs_marg,'Replace',false);
    end
end

% build the decoders and do the predictions
for d = 1:length(D)
    disp(['building decoders with randomly chosen dPCs as predictors for dataset ', num2str(D(d))]);
    for m = min_nbr_marg:max_nbr_marg
        for i = 1:nbr_rnd_iter
            dPCA_sign{d,m,i} = build_dPC_decoder( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'dpcs',             rand_dpcs{m,i}, ...
                            'xval_yn',          xval_yn, ...
                            'fold_size',        60, ...
                            'smooth_spikes',    true, ...
                            'dec_offset',       true, ...
                            'poly_order',       2, ...
                            'return_preds',     false, ...
                            'plot_yn',          false );
        end
    end
end
disp('...');


% ---------------------------------------------------------------
% 2. Plot

% a) summary plot per dataset
for d = 1:length(D)
    figure('units','normalized','outerposition',[0 0 1 1]);
    for m = 1:length(marg_set)
        subplot(1,length(marg_set),m)
        if ~xval_yn
            bar(dPCA_fit{d,m}.all.R2,'facecolor',[.6 .6 .6]),box off
            ylabel('R2')
        else
            hold on
            errorbar(mean(dPCA_fit{d,m}.mfxval.R2,2),std(dPCA_fit{d,m}.mfxval.R2,0,2),...
                'marker','none','color','k','linewidth',2,'linestyle','none')
            bar(mean(dPCA_fit{d,m}.mfxval.R2(:,2:end),2)),box off
            if m == 1, ylabel('xval R2'); end
        end
        xlim([0 size(dPCA_fit{d,m}.all.R2,1)+1])
        set(gca,'TickDir','out','FontSize',12)
        title( marg_leg(marg_set{m}) );
        ylim([0 1]);xlim([0 length(dPCA_fit{d,m}.all.R2)+1])
        if m == 1
            if ~xval_yn, ylabel('R2'), else ylabel('mfxval R2'), end
        end
    end
end



% ---------------------------------------------------------------
% 3. Stats

% % 2-way ANOVA comparing across sets of marginalizations and muscles for
% % each dataset
% % data has to be a 2D in which each color is a marginalization and blocks
% % of X rows are the values for each muscle (if using mfxval, it will be the
% % R2/r for each fold)
% 
% for d = 1:length(D)
% 
%     data_anova      = zeros(numel(dPCA_fit{d,1}.mfxval.R2),length(marg_set));
%     nbr_folds       = size(dPCA_fit{d,1}.mfxval.R2,2);
%     for m = 1:length(marg_set)
%         data_anova(:,m) = reshape(dPCA_fit{d,m}.mfxval.R2,[],1);
%     end
% 
%     [p{d} tbl{d}]   = anova2(data_anova,nbr_folds);
% end