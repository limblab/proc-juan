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

% use only good EMGs? --- THIS NEEDS TO BE BUILT IN BUILD_DP_DECODER, in
% the mean time is used for plotting and for the stats
good_emgs_only = true;
for d = 1:length(D)
    if good_emgs_only
        emg_chs{d} = datasets{D(d)}.chosen_emgs;
    else
        emg_chs{d} = length(datasets{D(d)}.binned_data{1}.emgguide);
    end
end


% marginalization legends for plotting fits
marg_leg    = {'task','target','time','task/tgt int'};

% combinations of marginalizations that will be compared
marg_set    = {1, 2, 3, 4, [2 4], [1 3], [2 3], [1 2 4], [1 2 3], [1 3 4], [1 2 3 4]};
num_margs   = length(marg_set);


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
nbr_rnd_iter = 100;

% get minimum and maximum number of input marginalizations
max_nbr_marg = max(cellfun(@(x) numel(x), marg_set));
min_nbr_marg = min(cellfun(@(x) numel(x), marg_set));
total_nbr_margs = max_nbr_marg - min_nbr_marg + 1;


% define random combinations of dPCs to use as predictors --will be the
% same for all the datasets
for n = min_nbr_marg:max_nbr_marg
    for i = 1:nbr_rnd_iter
        rand_dpcs{n,i} = datasample(1:mani_dim,n*dpcs_marg,'Replace',false);
    end
end

% build the decoders and do the predictions
for d = 1:length(D)
    disp(['building decoders with randomly chosen dPCs as predictors for dataset ', num2str(D(d))]);
    for n = min_nbr_marg:max_nbr_marg
        for i = 1:nbr_rnd_iter
            dPCA_sign{d,n,i} = build_dPC_decoder( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'dpcs',             rand_dpcs{n,i}, ...
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

% get mean, SD of the decoder predictions with randomly chosen inputs
for d = 1:length(D)
    for n = 1:size(dPCA_sign,2)
        dsc     = dPCA_sign(d,n,:);
        dsc     = reshape(dsc,[1,size(dPCA_sign,3)]);
        dsm     = cell2mat(cellfun(@(x) x.mfxval.R2, dsc,'UniformOutput', false));
        rnd_dec{d,n}.mn = mean(dsm,2);
        rnd_dec{d,n}.sd = std(dsm,0,2);
    end
end


% a) summary plot per dataset
for d = 1:length(D)
    figure('units','normalized','outerposition',[0 0 1 1]);
    for m = 1:num_margs
        
        % Bar plot for the predictions
        subplot(1,length(marg_set)+total_nbr_margs,m)
        if ~xval_yn
            bar(dPCA_fit{d,m}.all.R2(emg_chs{d}),'facecolor',[.6 .6 .6]),box off
            ylabel('R2')
        else
            hold on
            errorbar(mean(dPCA_fit{d,m}.mfxval.R2(emg_chs{d},:),2),...
                std(dPCA_fit{d,m}.mfxval.R2(emg_chs{d},:),0,2),...
                'marker','none','color','k','linewidth',1,'linestyle','none')
            bar(mean(dPCA_fit{d,m}.mfxval.R2(emg_chs{d},:),2)),box off
            if m == 1, ylabel('xval R2'); end
        end
        %xlim([0 size(dPCA_fit{d,m}.all.R2(emg_chs{d}),1)+1])
        set(gca,'TickDir','out','FontSize',12)
        title( marg_leg(marg_set{m}) );
        ylim([0 1]);xlim([0 length(dPCA_fit{d,m}.all.R2(emg_chs{d},:))+1])
        if m == 1
            if ~xval_yn, ylabel('R2'), else ylabel('mfxval R2'), end
        end
    end
    % plot predicitons with randomly selected inputs
    for m = 1:total_nbr_margs
        subplot(1,length(marg_set)+total_nbr_margs,length(marg_set)+m)
        if ~xval_yn
            disp('')
        else
            hold on
            errorbar(rnd_dec{d,m}.mn(emg_chs{d}),rnd_dec{d,m}.sd(emg_chs{d}),...
                'marker','none','color','k','linewidth',1,'linestyle','none')
            bar(rnd_dec{d,m}.mn(emg_chs{d}),'FaceColor','r'),box off
        end
        set(gca,'TickDir','out','FontSize',12)
        title({'rand. marg'; num2str(m)});
        ylim([0 1]);xlim([0 length(dPCA_fit{d,m}.all.R2(emg_chs{d},:))+1])
    end
end



% ---------------------------------------------------------------
% 3. Stats

% significance threshold
P_sign      = 0.05;


% get number of muscles per session
num_emgs_p_dataset = zeros(1,length(D));

for d = 1:length(D)
   num_emgs_p_dataset(d) = length(dPCA_fit{d,1}.all.r);
end

% matrix to store the results of the comparisons 
for d = 1:length(D)
    P_matrix_fits{d} = nan(num_margs*(num_margs-1),num_emgs_p_dataset(d));
end

% THINK OF A BETTER WAY TO ORGANIZE THIS !!!!

% a) Paired t-test to compare decoder predictions for each muscle and all
% pairs of combinations of tasks
for d = 1:length(D)
    for m = 1:num_margs
        
        % see what marginalizations to compare
        other_margs = setdiff(1:num_margs,m);
        
        if ~xval_yn
            error('NOT IMPLEMENTED WITHOUT MFXVAL'); % ToDo
        else
            for o = 1:length(other_margs)
                [~,p] = ttest(dPCA_fit{d,m}.mfxval.R2(emg_chs{d},:)',...
                            dPCA_fit{d,other_margs(o)}.mfxval.R2(emg_chs{d},:)');
                % store PVal
                P_matrix_fits{d}((m-1)*num_margs+o,:) = p;
            end
        end
        
    end
end


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