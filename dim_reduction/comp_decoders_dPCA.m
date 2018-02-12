%
% Compare EMG fits using decoders that take different sets
% (marginalizations) of dPCs as inputs
%

% clearvars -except dPCA_data datasets


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------

% check that the parallel pool is running, otherwise start it
gcp;

% ---------------------------------------------------------------
% Some definitions

% crossvalidate the predicitons?
xval_yn     = true;

% manifold dimensionality
mani_dim    = 12;

% number of dPCs per marginalization to use
dpcs_marg   = 2;

% use s.e.m. rather than s.d. for plotting
use_sem_yn  = false;

% fold size (s)
fold_size   = 30;

% order polynomial static non-linearity
poly_order  = 3;

% marginalization legends for plotting fits -- THIS DEFINES DATA IN ALL THE FIGURES
marg_leg    = {'task','target','time','task/tgt int'};

% combinations of marginalizations that will be compared
% marg_set    = {1, 2, 3, 4, [2 4], [1 3], [2 3], [1 2 4], [1 2 3], [1 3 4], [1 2 3 4]};
marg_set    = {1, 2, 3, 4, [1 2 3 4]};

% datasets that will be analyzed
D           = [1:3 7:9]; % Jaco and Jango 1D / 2D datasets

% use only good EMGs? --- THIS NEEDS TO BE BUILT IN BUILD_DP_DECODER, in
% the mean time is used for plotting and for the stats
good_emgs_only = true;


% -------------------------------------------------------------------------
% Some calculations based on the definitions above

for d = 1:length(D)
    if good_emgs_only
        emg_chs{d} = datasets{D(d)}.chosen_emgs;
    else
        emg_chs{d} = 1:length(datasets{D(d)}.binned_data{1}.emgguide);
    end
end


num_margs   = length(marg_set);

% get minimum and maximum number of input marginalizations
max_nbr_marg = max(cellfun(@(x) numel(x), marg_set));
min_nbr_marg = min(cellfun(@(x) numel(x), marg_set));
total_nbr_margs = max_nbr_marg - min_nbr_marg + 1;


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 0. Do dPCA of the neural data, if the data are not avaialble

if exist('dPCA_results','var')
    dPCA_data = dPCA_results;
end
if ~exist('dPCA_data','var')
    for d = 1:length(D)
        dPCA_data{d} = call_dPCA( datasets{D(d)}.stdata, mani_dim, false );
    end
end


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 1. Build decoders with N dPCs for the same MARGINALIZATIONS (from one
% or several marginalizations) as predictors

% Build decoders taking the leading 'dpcs_marg' dPCs per marginalization as
% predictors

for d = 1:length(D)
    disp(['building decoders with chosen marginalizations for dataset ', num2str(D(d))]);
    for m = 1:length(marg_set)
        %if m == 2, plot_flag = true; else plot_flag = false; end
        dPCA_fit{d,m} = build_dPC_decoder( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'margs',            marg_set{m}, ...
                            'nbr_dpcs_p_marg',  dpcs_marg, ... 
                            'xval_yn',          xval_yn, ...
                            'fold_size',        fold_size, ...
                            'smooth_spikes',    false, ...
                            'dec_offset',       true, ...
                            'poly_order',       poly_order, ...
                            'return_preds',     true, ...
                            'dec_p_task',       false, ...
                            'tgts_to_excl',     {[-1 -5 1 -7],[-1 7 1 5]}, ...
                            'plot_yn',          false );
    end
end
disp('...');



%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% 2. Build decoder with all the dPCs


for d = 1:length(D)
    disp(['building decoders with all the dPCs for dataset ', num2str(D(d))]);
    all_dPCA_fit{d} = build_dPC_decoder( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'dpcs',             1:mani_dim, ...
                            'xval_yn',          xval_yn, ...
                            'fold_size',        fold_size, ...
                            'smooth_spikes',    true, ...
                            'dec_offset',       true, ...
                            'poly_order',       poly_order, ...
                            'return_preds',     false, ...
                            'dec_p_task',       false, ...
                            'tgts_to_excl',     {[-1 -5 1 -7],[-1 7 1 5]}, ...
                            'plot_yn',          false );
end
disp('...');


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% Stats and some other definitions for the analysis

% Some definitions

ds_p_monk       = {[1 2 3],[4 5 6]}; % this is with respect to the dPC results
ds_all_p_monk   = {[1 2 3],[7 8 9]}; % this is wrt to the datasets struct
monk            = {'JacoB','Jango'};


% PERC NEURAL VARIANCE EXPLAINED per marginalization (the number of
% components)

perc_neural_var_marg    = cell(1,length(D));
for d = 1:6, 
    
    % find the two leading dPCs that correspond to each marginalization
    [~, marg_max ]      = max(dPCA_data{d}.expl_var.margVar);
    
    % calculate percentage variance explained
    for m = 1:length(dPCA_data{d}.expl_var.totalMarginalizedVar)
        pos_margs       = find( marg_max==m, dpcs_marg );
        perc_neural_var_marg{d}(m) = sum(dPCA_data{d}.expl_var.margVar(m,pos_margs))/100;
    end
    
%     % the old version of the code, the total variance corresponding to that
%     % marginalization
%     perc_neural_var_marg_2{d} = dPCA_data{d}.expl_var.totalMarginalizedVar/dPCA_data{d}.expl_var.totalVar; 
end



% ----------------------------
% POOL DATA OVER ALL MUSCLES FOR EACH SESSION

% R2 EMGs for all the muscles for each marginalization and session
R2_marg_all_musc        = cell(length(D),length(marg_leg));

for d = 1:length(D)
    for m = 1:length(marg_leg)
        R2_marg_all_musc{d,m} = reshape(dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:),1,[]);
    end
end


% R2 EMGs for all the muscles using all dPCs as predictors
R2_marg_all_musc_all_dPCs   = cell(length(D),1);

for d = 1:length(D)
    R2_marg_all_musc_all_dPCs{d} = reshape(all_dPCA_fit{d}.mfxval{1}.R2(emg_chs{d},:),1,[]);
end



% ----------------------------
% POOL DATA OVER ALL THE SESSIONS FOR EACH MONKEY SEPARATELY

% R2 EMGs POOLED OVER MONKEYS for all the sessions and muscles per marginalization
R2_marg_all_musc_per_monk   = cell(length(monk),length(marg_leg));

for k = 1:length(monk)
    % do per marginalization
    for m = 1:length(marg_leg)
        % pool all the session for each monkey
        for d = 1:numel(ds_p_monk{k})
            R2_marg_all_musc_per_monk{k,m} = [R2_marg_all_musc_per_monk{k,m}, ...
                                                R2_marg_all_musc{ds_p_monk{k}(d),m}];
        end
    end
end


% R2 EMGs POOLED OVER MONKEYS for all the sessions and muscle using all dPCs as predictors
R2_marg_all_musc_all_dPCs_per_monk  = cell(length(monk),1);

for k = 1:length(monk)
    for d = 1:numel(ds_p_monk{k})
        R2_marg_all_musc_all_dPCs_per_monk{k} = [R2_marg_all_musc_all_dPCs_per_monk{k}, ...
                                                R2_marg_all_musc_all_dPCs{ds_p_monk{k}(d)}];
    end
end


% ----------------------------
% POOL DATA OVER ALL THE SESSIONS AND MONKEYS


% R2 for each marginalization across all muscles and monkeys
R2_marg_all             = cell(1,length(marg_leg));

for m = 1:length(marg_leg)
    R2_marg_all{m}  = cell2mat(R2_marg_all_musc_per_monk(:,m)');                 
end


% R2 for across all muscles and monkeys using all dPCs as predictors
R2_marg_all_all_dPCs    = cell2mat(R2_marg_all_musc_all_dPCs');


% --------------------------------------------------------
% --------------------------------------------------------
% NORMALIZE THE DATA TO THE PREDICTIONS OBTAINED WITH ALL DPCS

norm_R2_marg_all_musc   = cell(length(D),length(marg_leg));

for d = 1:length(D)
    for m = 1:length(marg_leg)
        norm_R2_marg_all_musc{d,m} = reshape(dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:),1,[])./...
            reshape(all_dPCA_fit{d}.mfxval{1}.R2(emg_chs{d},:),1,[]);
    end
end



% ----------------------------
% POOL NORMALIZED DATA OVER ALL THE SESSIONS FOR EACH MONKEY SEPARATELY

% Norm R2 EMGs POOLED OVER MONKEYS for all the sessions and muscles per marginalization
norm_R2_marg_all_musc_per_monk = cell(length(monk),length(marg_leg));

for k = 1:length(monk)
    % do per marginalization
    for m = 1:length(marg_leg)
        % pool all the session for each monkey
        for d = 1:numel(ds_p_monk{k})
            norm_R2_marg_all_musc_per_monk{k,m} = [norm_R2_marg_all_musc_per_monk{k,m}, ...
                                                norm_R2_marg_all_musc{ds_p_monk{k}(d),m}];
        end
    end
end



% ----------------------------
% POOL NORMALIZED DATA OVER ALL THE SESSIONS AND MONKEYS


% Norm R2 for each marginalization across all muscles and monkeys
% R2 for each marginalization across all muscles and monkeys
norm_R2_marg_all            = cell(1,length(marg_leg));

for m = 1:length(marg_leg)
    norm_R2_marg_all{m}  = cell2mat(norm_R2_marg_all_musc_per_monk(:,m)');                 
end




% --------------------------------------------------------
% --------------------------------------------------------
% STATS PER SESSION


% for each marginalziation across all muscles

mn_R2_marg_all_musc     = zeros(length(D),length(marg_leg));
sd_R2_marg_all_musc     = zeros(length(D),length(marg_leg));
sem_R2_marg_all_musc    = zeros(length(D),length(marg_leg));

mn_norm_R2_marg_all_musc = zeros(length(D),length(marg_leg));
sd_norm_R2_marg_all_musc = zeros(length(D),length(marg_leg));
sem_norm_R2_marg_all_musc = zeros(length(D),length(marg_leg));



for d = 1:length(D)
    
    for m = 1:length(marg_leg)
        mn_R2_marg_all_musc(d,m)    = mean(R2_marg_all_musc{d,m});
        sd_R2_marg_all_musc(d,m)    = std(R2_marg_all_musc{d,m});
        sem_R2_marg_all_musc(d,m)   = sd_R2_marg_all_musc(d,m)/sqrt(length(R2_marg_all_musc{d,m}));
        
        mn_norm_R2_marg_all_musc(d,m) = mean(norm_R2_marg_all_musc{d,m});
        sd_norm_R2_marg_all_musc(d,m) = std(norm_R2_marg_all_musc{d,m});
        sem_norm_R2_marg_all_musc(d,m) = sd_norm_R2_marg_all_musc(d,m)/sqrt(length(norm_R2_marg_all_musc{d,m}));
    end
end


% using all dPCs across all muscles

mn_R2_marg_all_musc_all_dPCs    = cellfun(@(x) mean(x), R2_marg_all_musc_all_dPCs );
sd_R2_marg_all_musc_all_dPCs    = cellfun(@(x) std(x), R2_marg_all_musc_all_dPCs );
sem_R2_marg_all_musc_all_dPCs   = cellfun(@(x) std(x)/sqrt(length(x)), R2_marg_all_musc_all_dPCs );



    
% ----------------------------
% STATS PER MONKEY
    

% Neural variance explained per monkey
for k = 1:numel(monk)
    neural_marg_mn{k} = mean( cell2mat(perc_neural_var_marg(ds_p_monk{k})') );
    neural_marg_sd{k} = std( cell2mat(perc_neural_var_marg(ds_p_monk{k})') );
    neural_marg_sem{k} = std( cell2mat(perc_neural_var_marg(ds_p_monk{k})') )/...
                            sqrt(size(cell2mat(perc_neural_var_marg(ds_p_monk{k})'),1));
end


% R2 for each marginalziation across all muscles for each monkey

mn_R2_marg_all_musc_per_monk    = cellfun(@(x) mean(x), R2_marg_all_musc_per_monk(:,:) );
sd_R2_marg_all_musc_per_monk    = cellfun(@(x) std(x), R2_marg_all_musc_per_monk(:,:) );
sem_R2_marg_all_musc_per_monk   = cellfun(@(x) std(x)/sqrt(numel(x)), R2_marg_all_musc_per_monk(:,:) );


% R2 using all dPCs across all muscles for each monkey

mn_R2_marg_all_musc_all_dPCs_per_monk   = cellfun( @(x) mean(x), R2_marg_all_musc_all_dPCs_per_monk);
sd_R2_marg_all_musc_all_dPCs_per_monk   = cellfun( @(x) std(x), R2_marg_all_musc_all_dPCs_per_monk);
sem_R2_marg_all_musc_all_dPCs_per_monk  = cellfun( @(x) std(x)/sqrt(numel(x)), R2_marg_all_musc_all_dPCs_per_monk);


% Norm. R2 for each marginalziation across all muscles for each monkey

mn_norm_R2_marg_all_musc_per_monk       = cellfun(@(x) mean(x), norm_R2_marg_all_musc_per_monk(:,:) );
sd_norm_R2_marg_all_musc_per_monk       = cellfun(@(x) std(x), norm_R2_marg_all_musc_per_monk(:,:) );
sem_norm_R2_marg_all_musc_per_monk      = cellfun(@(x) std(x)/sqrt(numel(x)), norm_R2_marg_all_musc_per_monk(:,:) );


% ----------------------------
% STATS FOR BOTH MONKEYS TOGETHER

mn_R2_marg_all          = cellfun( @(x) mean(x), R2_marg_all );
sd_R2_marg_all          = cellfun( @(x) std(x), R2_marg_all );
sem_R2_marg_all         = cellfun( @(x) std(x)/sqrt(numel(x)), R2_marg_all );

mn_R2_marg_all_all_dPCs = mean(R2_marg_all_all_dPCs);
sd_R2_marg_all_all_dPCs = std(R2_marg_all_all_dPCs);
sem_R2_marg_all_all_dPCs = std(R2_marg_all_all_dPCs)/sqrt(numel(R2_marg_all_all_dPCs));

mn_norm_R2_marg_all     = cellfun( @(x) mean(x), norm_R2_marg_all );
sd_norm_R2_marg_all     = cellfun( @(x) std(x), norm_R2_marg_all );
sem_norm_R2_marg_all    = cellfun( @(x) std(x)/sqrt(numel(x)), norm_R2_marg_all );



% colors for the plots

col_all_dpcs    = [.9 .5 .9];
marg_cols       = [.4 .6 .5; .7 .1 .1; .5 .65 .9; 1 .6 .3];



%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 2. Plot

% a) summary plot per dataset: each muscle for each marginalization, and
% all the dPCs

for d = 1:length(D)
    
    figure
    
    for m = 1:num_margs    
        
        % Bar plot for the predictions
        subplot(1,length(marg_set),m)
        if ~xval_yn
            bar(dPCA_fit{d,m}.all.R2(emg_chs{d}),'facecolor',[.6 .6 .6]),box off
            ylabel('R2')
        else
            hold on
            eb      = std(dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:),0,2);
            if use_sem_yn
                eb = eb/sqrt(size(dPCA_fit{d,m}.mfxval{1}.R2,2));
            end
            errorbar(mean(dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:),2),...
                eb,...
                'marker','none','color','k','linewidth',1,'linestyle','none')
            tb = bar(mean(dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:),2)); box off
            if m <= size(marg_cols,1)
                set(tb,'FaceColor',marg_cols(m,:),'EdgeColor',marg_cols(m,:));
            elseif m == size(marg_cols,1)+1
                % set(tb,'FaceColor',[155 110 175]/255,'EdgeColor',[155 110 175]/255);
                set(tb,'FaceColor',col_all_dpcs,'EdgeColor',col_all_dpcs);
            end
            if m == 1, ylabel('xval R2'); end
        end
        
        %xlim([0 size(dPCA_fit{d,m}.all.R2(emg_chs{d}),1)+1])
        set(gca,'TickDir','out','FontSize',12)
        title( marg_leg(marg_set{m}) );
        ylim([0 1]);xlim([0 length(dPCA_fit{d,m}.all{1}.R2(emg_chs{d},:))+1])
        if m == 1
            if ~xval_yn, ylabel('R2'), else ylabel('mfxval R2'), end
        end
    end
end



%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% PLOT ALL MUSCLES TOGETHER for each session, 
% for each of the 4 DPCA MARGINALIZATIONS, and all dPCs

nbr_subplots    = num_margs;


for d = 1:length(D)
    figure, hold on
    for s = 1:nbr_subplots
        if s <= 4
            % Bar plot with the dPCs
            % do one single bar for all the muscles
            if use_sem_yn, eb = sem_R2_marg_all_musc(d,s); else eb = sd_R2_marg_all_musc(d,s); end
            errorbar(s,mn_R2_marg_all_musc(d,s), eb,...
                'marker','none','color',marg_cols(s,:),'linewidth',1,'linestyle','none')
            bar(s,mn_R2_marg_all_musc(d,s),'FaceColor',marg_cols(s,:)),box off
        end
        % Predictions with all the dPCs
        if  s == 5
            % do one single bar for all the muscles
            if use_sem_yn, eb = sem_R2_marg_all_musc_all_dPCs(d); else eb = sd_R2_marg_all_musc_all_dPCs(d); end
            errorbar(s,mn_R2_marg_all_musc_all_dPCs(d), eb,...
                'marker','none','color',col_all_dpcs,'linewidth',1,'linestyle','none')
            bar(s,mn_R2_marg_all_musc_all_dPCs(d),'FaceColor',col_all_dpcs),box off
        end
    end
    ylim([0 1]);xlim([0 nbr_subplots+1])
    ylabel('R^2 EMG Predictions')
    set(gca,'TickDir','out','FontSize',12)
    set(gca,'XTick',1:5,'XTickLabel',{'task','target','dynamics','task/target','all dPCs'},...
        'XTickLabelRotation',45)
end


%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT ALL MUSCLES TOGETHER for the 3 sessions by each monkey, 
% for each of the 4 DPCA MARGINALIZATIONS, and all dPCs


% 1 plot per monkey
for k = 1:numel(ds_p_monk)
    
    figure, hold on
    % do for each subplot
    for s = 1:5
        
        % for the dPC marginalizations
        if s <= 4

            if use_sem_yn, eb = sem_R2_marg_all_musc_per_monk(k,s); else, eb = sd_R2_marg_all_musc_per_monk(k,s); end
            % and plot
            errorbar(s,mn_R2_marg_all_musc_per_monk(k,s), eb,...
                'marker','none','color',marg_cols(s,:),'linewidth',1,'linestyle','none')
            bar(s,mn_R2_marg_all_musc_per_monk(k,s),'FaceColor',marg_cols(s,:)),box off
        end
        
         % Predictions with all the dPCs
        if s == 5
            
            if use_sem_yn, eb = sem_R2_marg_all_musc_all_dPCs_per_monk(k); else, eb = sd_R2_marg_all_musc_all_dPCs_per_monk(k); end
            % and plot
            errorbar(s,mn_R2_marg_all_musc_all_dPCs_per_monk(k), eb,...
                'marker','none','color',col_all_dpcs,'linewidth',1,'linestyle','none')
            bar(s,mn_R2_marg_all_musc_all_dPCs_per_monk(k),'FaceColor',col_all_dpcs),box off
        end
       
    end
    ylim([0 1]);xlim([0 s+1])
    ylabel('R^2 EMG Predictions')
    set(gca,'TickDir','out','FontSize',12)
    set(gca,'XTick',1:5,'XTickLabel',{'task','target','dynamics','task/target','all dPCs'},'XTickLabelRotation',45)
    title(monk{k})
end
 


%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% SAME AS PREVIOUS PLOT FOR BOTH MONKEYS TOGETHER


figure, hold on
for s = 1:nbr_subplots
	% for the dPC marginalizations
    if s <= 4
        
        if use_sem_yn, eb = sem_R2_marg_all(s); else eb = sd_R2_marg_all(s); end
        % and plot
        errorbar(s,mn_R2_marg_all(s), eb,...
                'marker','none','color',marg_cols(s,:),...
                'linewidth',1,'linestyle','none')
        bar(s,mn_R2_marg_all(s),'FaceColor',marg_cols(s,:))
    end
    
    % all dPCs
     if s == 5
         
        if use_sem_yn, eb = sem_R2_marg_all_all_dPCs; else eb = sd_R2_marg_all_all_dPCs; end
        % and plot
        errorbar(s,mn_R2_marg_all_all_dPCs, eb,...
                'marker','none','color',col_all_dpcs,...
                'linewidth',1,'linestyle','none')
        bar(s,mn_R2_marg_all_all_dPCs,'FaceColor',col_all_dpcs)
     end
end
ylim([0 1]); xlim([0 nbr_subplots+1]); box off
ylabel('R^2 EMG Predictions')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:nbr_subplots,'XTickLabel',{'task','target','dynamics','task/target','all dPCs'},'XTickLabelRotation',45)




%% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% NORMALIZED PLOTS


%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% a) summary plot per dataset: each muscle for each marginalization, and
% all the dPCs


% do not plot the marginalization that corresponds to all the co-variates
% (all_dPCs)
if sum( cellfun(@(x) numel(x), marg_set) == 4 ) > 0
    num_subplots = length(marg_set)-1;
else
    num_subplots = length(marg_set);
end


for d = 1:length(D)
    
    figure
    
    for m = 1:num_subplots   
        
        % Bar plot for the predictions
        subplot(1,num_subplots,m)
        if ~xval_yn
            bar(dPCA_fit{d,m}.all.R2(emg_chs{d}),'facecolor',[.6 .6 .6]),box off
            ylabel('R2')
        else
            hold on
            data    = dPCA_fit{d,m}.mfxval{1}.R2(emg_chs{d},:)./...
                        all_dPCA_fit{d}.mfxval{1}.R2(emg_chs{d},:);
            mn_data = mean(data,2);
            eb      = std(data,0,2);
            if use_sem_yn, eb  = eb/sqrt(size(data,2)); end
            errorbar( mn_data, eb, 'marker','none','color','k',...
                    'linewidth',1,'linestyle','none' )
            tb      = bar(mn_data); box off
            set(tb,'FaceColor',marg_cols(m,:),'EdgeColor',marg_cols(m,:));

            if m == 1, ylabel('xval Norm. R2'); end
        end
        
        %xlim([0 size(dPCA_fit{d,m}.all.R2(emg_chs{d}),1)+1])
        set(gca,'TickDir','out','FontSize',12)
        title( marg_leg(marg_set{m}) );
        ylim([0 1]);xlim([0 length(dPCA_fit{d,m}.all{1}.R2(emg_chs{d},:))+1])
        if m == 1
            if ~xval_yn, ylabel('R2'), else ylabel('mfxval Norm. R2'), end
        end
    end
end

clear eb




%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% PLOT ALL MUSCLES TOGETHER for each session, 
% for each of the 4 DPCA MARGINALIZATIONS, and all dPCs


% do not plot the marginalization that corresponds to all the co-variates
% (all_dPCs)
if sum( cellfun(@(x) numel(x), marg_set) == 4 ) > 0
    num_subplots = length(marg_set)-1;
else
    num_subplots = length(marg_set);
end


for d = 1:length(D)
    figure, hold on
    for s = 1:num_subplots
        if s <= 4
            % Bar plot with the dPCs
            % do one single bar for all the muscles
            if use_sem_yn, eb = sem_norm_R2_marg_all_musc(d,s); else eb = sd_norm_R2_marg_all_musc(d,s); end
            errorbar(s,mn_norm_R2_marg_all_musc(d,s), eb,...
                'marker','none','color',marg_cols(s,:),'linewidth',1,'linestyle','none')
            bar(s,mn_norm_R2_marg_all_musc(d,s),'FaceColor',marg_cols(s,:)),box off
        end

    end
    ylim([0 1]);xlim([0 num_subplots+1])
    ylabel('Norm. R^2 EMG Predictions')
    set(gca,'TickDir','out','FontSize',12)
    set(gca,'XTick',1:num_subplots,'XTickLabel',{'task','target','dynamics','task/target'},...
        'XTickLabelRotation',45)
end



%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT ALL MUSCLES TOGETHER for the 3 sessions by each monkey, 
% for each of the 4 DPCA MARGINALIZATIONS


% 1 plot per monkey
for k = 1:numel(ds_p_monk)
    
    figure, hold on
    % do for each subplot
    for s = 1:num_subplots
        
        % for the dPC marginalizations
        if use_sem_yn, eb = sem_norm_R2_marg_all_musc_per_monk(k,s); else eb = sd_norm_R2_marg_all_musc_per_monk(k,s); end
        % and plot
        errorbar(s,mn_norm_R2_marg_all_musc_per_monk(k,s), eb,...
            'marker','none','color',marg_cols(s,:),'linewidth',1,'linestyle','none')
        bar(s,mn_norm_R2_marg_all_musc_per_monk(k,s),'FaceColor',marg_cols(s,:)),box off
       
    end
    ylim([0 1]);xlim([0 num_subplots+1])
    ylabel('Norm. R^2 EMG Predictions')
    set(gca,'TickDir','out','FontSize',12)
    set(gca,'XTick',1:num_subplots,'XTickLabel',{'task','target','dynamics','task/target'},'XTickLabelRotation',45)
    title(monk{k})
end
 


%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% SAME AS PREVIOUS PLOT FOR BOTH MONKEYS TOGETHER


figure, hold on
for s = 1:num_subplots
	% for the dPC marginalizations
    if s <= 4
        
        if use_sem_yn, eb = sem_norm_R2_marg_all(s); else eb = sd_norm_R2_marg_all(s); end
        % and plot
        errorbar(s,mn_norm_R2_marg_all(s), eb,...
                'marker','none','color',marg_cols(s,:),...
                'linewidth',1,'linestyle','none')
        bar(s,mn_norm_R2_marg_all(s),'FaceColor',marg_cols(s,:))
    end
end
ylim([0 1]); xlim([0 num_subplots+1]); box off
ylabel('Norm. R^2 EMG Predictions')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:num_subplots,'XTickLabel',{'task','target','dynamics','task/target'},'XTickLabelRotation',45)



%% -----------------------------------------------------------------------
% ---------------------------------------------------------------
% SAME AS PREVIOUS PLOT BUT ALSO INCLUDING PREDICTIONS FROM ALL 12 DPCS,
% AND NOT NORMALIZED BUT "RAW" R^2


figure, hold on
for s = 1:num_subplots
	% for the dPC marginalizations
    if s <= 4
        
        if use_sem_yn, eb = sem_R2_marg_all(s); else eb = sd_R2_marg_all(s); end
        % and plot
        errorbar(s,mn_R2_marg_all(s), eb,...
                'marker','none','color',marg_cols(s,:),...
                'linewidth',1,'linestyle','none')
        bar(s,mn_R2_marg_all(s),'FaceColor',marg_cols(s,:))
    end
end
% Add results all the dPCs
if use_sem_yn, eb = sem_R2_marg_all_all_dPCs; else eb = sd_R2_marg_all_all_dPCs; end
errorbar(num_subplots+1,mn_R2_marg_all_all_dPCs,eb,...
    'marker','none','color',col_all_dpcs,...
	'linewidth',1,'linestyle','none')
bar(num_subplots+1,mn_R2_marg_all_all_dPCs,'FaceColor',col_all_dpcs);
ylim([0 1]); xlim([0 num_subplots+2]); box off
ylabel('Norm. R^2 EMG Predictions')
set(gca,'TickDir','out','FontSize',12)
set(gca,'XTick',1:num_subplots+1,'XTickLabel',{'task','target','dynamics','task/target','all dPCs'},'XTickLabelRotation',45)




%% -----------------------------------------------------------------------
% ------------------------------------
% PLOT EXPLAINED NEURAL VARIANCE VS. EMG PREDICTIONS


% For all 4 marginalizations of the dPCs

figure,hold on
for k = 1:numel(monk)
    for m = 1:4
        if k == 1 % Jaco
            plot(neural_marg_mn{k}(m),mn_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor','none')
        else
            plot(neural_marg_mn{k}(m),mn_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor',marg_cols(m,:))
        end
        plot([neural_marg_mn{k}(m)-neural_marg_sd{k}(m),...
                neural_marg_mn{k}(m)+neural_marg_sd{k}(m)],...
                [mn_R2_marg_all_musc_per_monk(k,m), mn_R2_marg_all_musc_per_monk(k,m)],...
                'color',marg_cols(m,:))
        plot([neural_marg_mn{k}(m), neural_marg_mn{k}(m)],...
            [mn_R2_marg_all_musc_per_monk(k,m)-sd_R2_marg_all_musc_per_monk(k,m),...
            mn_R2_marg_all_musc_per_monk(k,m)+sd_R2_marg_all_musc_per_monk(k,m)],...
            'color',marg_cols(m,:))
    end
end


% % For all the dPCs
% if exist('all_dPCA_fit','var')
%     for k = 1:numel(monk)
%         if k == 1
%               plot(1,mean(all_dpcs{k}),'s','markersize',14,...
%                     'linewidth',1.5,'color',col_all_dpcs,'MarkerFaceColor','none')
%         else
%               plot(1,mean(all_dpcs{k}),'s','markersize',14,...
%                     'linewidth',1.5,'color',col_all_dpcs,'MarkerFaceColor',col_all_dpcs)
%         end
%             if ~use_sem_yn
%                 plot([1 1],...
%                     [mean(all_dpcs{k})+std(all_dpcs{k}),...
%                     mean(all_dpcs{k})-std(all_dpcs{k})],...
%                     'color',col_all_dpcs)
%             else
%                 plot([1 1],...
%                     [mean(all_dpcs{k})+std(all_dpcs{k})/sqrt(length(all_dpcs{k})),...
%                     mean(all_dpcs{k})-std(all_dpcs{k})/sqrt(length(all_dpcs{k}))],...
%                     'color',col_all_dpcs)
%             end
%     end
% end
ylim([0 .6]);xlim([0 .4])
ylabel('R^2 EMG')
xlabel('Neural varance expl. (%)')
set(gca,'TickDir','out','FontSize',12)
   


%% -----------------------------------------------------------------------
% ------------------------------------
% SAME AS PREVIOUS PLOT BUT WITH NORMALIZED EMG PREDICTIONS

figure,hold on
for k = 1:numel(monk)
    for m = 1:4
        if k == 1 % Jaco
            plot(neural_marg_mn{k}(m),mn_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor','none')
        else
            plot(neural_marg_mn{k}(m),mn_norm_R2_marg_all_musc_per_monk(k,m),'s','markersize',14,...
                'linewidth',1.5,'color',marg_cols(m,:),'MarkerFaceColor',marg_cols(m,:))
        end
        plot([neural_marg_mn{k}(m)-neural_marg_sd{k}(m),...
                neural_marg_mn{k}(m)+neural_marg_sd{k}(m)],...
                [mn_norm_R2_marg_all_musc_per_monk(k,m), mn_norm_R2_marg_all_musc_per_monk(k,m)],...
                'color',marg_cols(m,:))
        plot([neural_marg_mn{k}(m), neural_marg_mn{k}(m)],...
            [mn_norm_R2_marg_all_musc_per_monk(k,m)-sd_norm_R2_marg_all_musc_per_monk(k,m),...
            mn_norm_R2_marg_all_musc_per_monk(k,m)+sd_norm_R2_marg_all_musc_per_monk(k,m)],...
            'color',marg_cols(m,:))
    end
end


% % For all the dPCs
% if exist('all_dPCA_fit','var')
%     for k = 1:numel(monk)
%         if k == 1
%               plot(1,mean(all_dpcs{k}),'s','markersize',14,...
%                     'linewidth',1.5,'color',col_all_dpcs,'MarkerFaceColor','none')
%         else
%               plot(1,mean(all_dpcs{k}),'s','markersize',14,...
%                     'linewidth',1.5,'color',col_all_dpcs,'MarkerFaceColor',col_all_dpcs)
%         end
%             if ~use_sem_yn
%                 plot([1 1],...
%                     [mean(all_dpcs{k})+std(all_dpcs{k}),...
%                     mean(all_dpcs{k})-std(all_dpcs{k})],...
%                     'color',col_all_dpcs)
%             else
%                 plot([1 1],...
%                     [mean(all_dpcs{k})+std(all_dpcs{k})/sqrt(length(all_dpcs{k})),...
%                     mean(all_dpcs{k})-std(all_dpcs{k})/sqrt(length(all_dpcs{k}))],...
%                     'color',col_all_dpcs)
%             end
%     end
% end
ylim([0 1]);xlim([0 .5])
ylabel('Norm. R^2 EMG')
xlabel('Neural varance expl. (%)')
set(gca,'TickDir','out','FontSize',12)
   