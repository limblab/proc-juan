%
%
%
% ZSCORING DOESN'T WORK WELL


function results = classify_across_days( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
in                  = [];
out                 = [];
method              = 'Bayes'; % 'NN' 'Bayes'
% win                     = [];
hist_bins           = [];
win_size            = [];
n_folds             = 5; % 'cca' or 'procrustes'
manifold            = [];
mani_dims           = 1:6;
zsc                 = false;


if nargin > 1, assignParams(who,params); end % overwrite defaults



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(in), error('Must provide decoder inputs'); end
if isempty(out), error('Must provide decoder output'); end
if isempty(manifold), error('Must provide decoder output'); end
if hist_bins<0, error('Must provide decoder output'); end
% ToDo: Implement some more checks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STUFF

% get all pairs of sessions
sessions            = unique({td.date});
n_sessions          = length(sessions);
comb_sessions       = nchoosek(1:n_sessions,2);
n_comb_sessions     = size(comb_sessions,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAYING WITH THINGS -- NEED TO BE MOVED TO PARAMETERS WHEN THIS IS A REAL
% FUNCTION
td_orig             = td;

% downsample again
%win_size                = 0.150; %0.150;
down_rate           = round(win_size/td_orig(1).bin_size);

% Only keep the last 'win_size' s of the trial
td                  = trimTD(td_orig,{'end',-down_rate},{'end',0});

% Make one big bin (a la Santhanam et al JNP 2007)
td                	= binTD(td,down_rate);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

disp(['Testing classifier stability using ' in ' inputs']);
disp(['Type: ' method]);
disp(['Classifying: ' out]);



% For all pairs of sessions, build a model on session 1 and test it on
% session 2
for c = 1:n_comb_sessions


    % =====================================================================
    % CLASSIFIER BASED ON ALIGNED LATENT ACTIVITY
    % =====================================================================
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALIGN LATENT ACTIVITY -- HAS TO BE DONE WITH THE ORIGINAL DATA,
    %                   BECAUSE THAT'S WHAT WE USE FOR THE OTHER ANALYSES
    
    % a) get two TD structs, one per session to compare
    trials1         = getTDidx(td_orig,'date',sessions{comb_sessions(c,1)});
    trials2         = getTDidx(td_orig,'date',sessions{comb_sessions(c,2)});


    % b) "align" dynamics with CCA
    cca_info(c)     = compDynamics( td_orig, manifold, trials1, trials2, mani_dims );


    % c) Add latent activity to the tdi structs
    td1             = td(trials1);
    td2             = td(trials2);
    
    for t = 1:length(trials1)
        proj1       = td(trials1(t)).PMd_pca(:,mani_dims) * cca_info(c).A;
        proj2       = td(trials2(t)).PMd_pca(:,mani_dims) * cca_info(c).B;
        td1(t).([manifold '_align']) = proj1;
        td2(t).([manifold '_align']) = proj2;
    end
    

    % d) Duplicate and shift to have history
    if hist_bins > 0
        td1      	= dupeAndShift(td1,{[manifold '_align'],hist_bins});
        td2        	= dupeAndShift(td2,{[manifold '_align'],hist_bins});
    end


    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPARE SOME EXTRA THINGS

    % overwrite classifier input
    if hist_bins > 0
        clas_input  = [manifold '_align_shift'];
    else
        clas_input 	= [manifold '_align'];
    end
    
    
    % TEMP: Do for each bin (after dupe and shifting)
    n_clas_vars     = size(td1(1).(clas_input),2);


    % TEMP: do not need this either
    b               = 1;
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CROSS-VALIDATED WITHIN-DAY PREDICTIONS FOR DAY 1
    %       - This is currently done without optimizing the classifier
    
    % a) To cross-validate, shuffle the trials because they are sorted by
    % target location
    ishuf           = randperm(length(td1));
    bins_p_fold     = floor(length(td1)/n_folds);

    % To store the results
    perf_train      = zeros(1,n_folds);
    perf_test       = zeros(1,n_folds);
    
    % b) Do
    for f = 1:n_folds
        
        % i) Retrieve the training and testing bins
        itest       = ishuf( (f-1)*bins_p_fold + 1 : f*bins_p_fold );
        itrain      = ishuf(~ismember(ishuf,itest));

        % ii) Prepare the data for the classifier -- Dimensions: 
        % X: trial x datapoints; Y: trial x target
        Xtrain      = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                        td1(itrain)', 'uniformoutput', false ) );
        Ytrain      = [td1(itrain).target_direction];
        
        if zsc
            Xtrain  = zscore(Xtrain);
        end
        
        % iii) Train the classifier
        switch method
            case 'NN'
                clas_xval = fitcknn(Xtrain,Ytrain);
            case 'Bayes'
                clas_xval = fitcnb(Xtrain,Ytrain);
            otherwise
                error('Only NN implemented so far');
        end

        % iv) Test performance on the testing set
        Xtest       = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td1(itest)', 'uniformoutput', false ) );            
        Ytest       = [td1(itest).target_direction];

        if zsc
            Xtest   = zscore(Xtest);
        end
            
        pred_test   = predict(clas_xval,Xtest)';
        
        % Performance (% succesfully classified targets)
        perf_test1(f) = ( length(pred_test) - sum( ( pred_test - Ytest ~= 0 ) ) ) / length(pred_test)*100;
        
        err_test1(f) = mean(abs(Ytest - pred_test));
        

        % v) For reference, compute performance on the training set
        pred_train  = predict(clas_xval,Xtrain)';
        
        perf_train(f) = ( length(pred_train) - sum( ( pred_train - Ytrain ~= 0 ) ) ) / length(pred_train)*100;
    end
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRAIN CLASSIFIER ON DAY 1

    
    % a) Prepare the data
    X1              = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                        td1', 'uniformoutput', false ) );
    Y1              = [td1.target_direction];    

    if zsc % z-score recommended for Bayes
        X1          = zscore(X1);
    end
    
    
    % b) Train classifier
    switch method
        case 'NN'
            clas = fitcknn(X1,Y1);
        case 'Bayes'
%             clas = fitcnb(X1,Y1);            
            clas = fitcnb(X1,Y1,'DistributionNames','kernel', ...
                'Kernel','triangle', ...
                'Support','unbounded', ...
                'OptimizeHyperparameters',{'Width'},...
                'HyperparameterOptimizationOptions',struct('ShowPlots',false));
        otherwise
            error('Only NN and Bayes implemented so far');
    end
    
    
    % c) Test performance on training set -- Remember that the optimal
    %       parameters where identified with cross-validation, thus these
    %       predictions are cross-validated
    pred_within     = predict(clas,X1)';
    perf_within     = ( length(pred_within) - sum( ( pred_within - Y1 ~= 0 ) ) ) / length(pred_within)*100;
    
    
    % d) Compute classif error (degrees)
    err_within      = mean(abs(Y1-pred_within));
    
    

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEST ON DAY 2
    
    
    % a) Prepare data from day 2
    X2              = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td2', 'uniformoutput', false ) );            
    Y2              = [td2.target_direction];

    if zsc % z-score recommended for Bayes
        X2          = zscore(X2);
    end
    

    % b) Test performance on day 2
    pred_across     = predict(clas,X2)';
    perf_across     = ( length(pred_across) - sum( (pred_across - Y2 ~= 0) ) ) / length(pred_across)*100;
   
    
    % c) Classifier error
    err_across      = mean(abs(Y2-pred_across));

    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WITHIN-DAY DECODER FOR DAY 2, TO NORMALIZE THE TEST RESULTS

    % a) Prepare the data
    
    clas_input_within2 = clas_input(1:end-6); % decode based on the raw latent activity
    
            
    % b) To cross-validate, shuffle the trials because they are sorted by
    % target location
    ishuf2          = randperm(length(td2));
    bins_p_fold2    = floor(length(td2)/n_folds);

    % To store the results
    perf_train2     = zeros(1,n_folds);
    perf_test2      = zeros(1,n_folds);
    
    
    % c) Do
    for f = 1:n_folds
        
        % i) Retrieve the training and testing bins
        itest2      = ishuf2( (f-1)*bins_p_fold2 + 1 : f*bins_p_fold2 );
        itrain2     = ishuf2(~ismember(ishuf2,itest2));

        % ii) Prepare the data for the classifier -- Dimensions: 
        % X: trial x datapoints; Y: trial x target
        Xtrain2     = cell2mat( arrayfun( @(x) x.(clas_input_within2)(b,1:n_clas_vars), ...
                        td2(itrain2)', 'uniformoutput', false ) );
        Ytrain2     = [td2(itrain2).target_direction];
        
        if zsc
            Xtrain2 = zscore(Xtrain2);
        end
        
        % iii) Train the classifier
        switch method
            case 'NN'
                clas_xval2 = fitcknn(Xtrain2,Ytrain2);
            case 'Bayes'
                clas_xval2 = fitcnb(Xtrain2,Ytrain2);
            otherwise
                error('Only NN implemented so far');
        end

        % iv) Test performance on the testing set
        Xtest2      = cell2mat( arrayfun( @(x) x.(clas_input_within2)(b,1:n_clas_vars), ...
                            td2(itest2)', 'uniformoutput', false ) );            
        Ytest2       = [td2(itest2).target_direction];

        if zsc
            Xtest2  = zscore(Xtest2);
        end
            
        pred_test2  = predict(clas_xval2,Xtest2)';
        
        % v) Performance (% succesfully classified targets)
        perf_test2(f) = ( length(pred_test2) - sum( ( pred_test2 - Ytest2 ~= 0 ) ) ) / length(pred_test2)*100;
        err_test2(f) = mean(abs(Ytest2-pred_test2));
        
        % For reference, compute performance on the training set
        pred_train2  = predict(clas_xval2,Xtrain2)';
        perf_train2(f) = ( length(pred_train2) - sum( ( pred_train2 - Ytrain2 ~= 0 ) ) ) / length(pred_train2)*100;
    end
    
    
%     % d) Test performance on training set 
%     pred_within2    = predict(clas_xval2,Xwithin2)';
%     perf_within2    = ( length(pred_within2) - sum( ( pred_within2 - Y2 ~= 0 ) ) ) / length(pred_within2)*100;
%     err_within2     = mean(abs(Y2-pred_within2));
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLASSIFIERS BASED ON SPIKES, WITHIN AND ACROSS DAY
    
    % ToDo
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save results
    
    
    % Optimized classifier, not cross-validated on day 1
    res.perf_within_opt_noxval1(c,b) = perf_within;
    res.err_within_opt_noxval1(c,b) = err_within;
    
    % Cross-validated on day 1
    res.perf_within_xval1(c,b,:) = perf_test1;
    res.err_within_xval1(c,b,:) = err_test1;
%     % not-cross-validated on day 1
%     res.perf_within_xval_on_train_set(c,b,:) = perf_train;
        
    % Cross-validated on day 2
    res.perf_within_xval2(c,b,:) = perf_test2;
    res.err_within_xval2(c,b,:) = err_test2;

    % Across day classifier
    res.perf_across(c,b) = perf_across;
    res.err_across(c,b) = err_across;
    
    % For the normalized predictions
    res.norm_perf_across(c,b) = perf_across/mean(perf_test2)*100;
    
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY STATS AND PLOTS


% Compute mean training performance per day
res.perf_within_xval1_m  = mean(res.perf_within_xval1,3);
res.perf_within_xval2_m  = mean(res.perf_within_xval2,3);

res.err_within_xval1_m  = mean(res.err_within_xval1,3);
res.err_within_xval2_m  = mean(res.err_within_xval2,3);


% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end
res.diff_days           = diff_days;
res.comb_sessions       = comb_sessions;


results                 = res;




% Histograms of metrics
hist_res                = 5;
x_hist                  = 0:hist_res:100+hist_res;
hist_within             = histcounts(res.perf_within_xval1_m(:,end),x_hist)/numel(res.perf_within_xval1_m)*100;
hist_across             = histcounts(res.perf_across(:,end),x_hist)/numel(res.perf_across)*100;

x_hist_norm             = 0:hist_res:100+hist_res+50;
hist_norm_across        = histcounts(res.norm_perf_across(:,end),x_hist_norm)/numel(res.norm_perf_across)*100;

Pchance                 = 100/numel(unique([td1.target_direction]));


% linear fits
x_fit                   = [0 max(diff_days)];
lfacross                = polyfit(diff_days,res.perf_across',1);
lfwithin2               = polyfit(diff_days,res.perf_within_xval2_m',1);
y_across                = polyval(lfacross,x_fit);
y_within2               = polyval(lfwithin2,x_fit);



% Fig

f1 = figure('units','normalized','outerposition',[0.11 0.1 0.9 0.9]);

subplot(331), hold on
errorbar(mean(res.perf_within_xval1_m,1),std(res.perf_within_xval1_m,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(res.perf_within_xval1_m,1),'FaceColor',[.7 .7 .7])
plot([0 size(res.perf_within_xval1_m,2)+1],[100 100], '-.', 'color', [.7 .7 .7], 'linewidth',1.5)
plot([0 size(res.perf_within_xval1_m,2)+1],[Pchance Pchance], '-.k', 'linewidth',1.5)
ylim([0 110])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Parameter b'),ylabel('Within-session accuracy (%)')
title(['Classifier: ' method])

subplot(332), hold on
errorbar(mean(res.perf_across,1),std(res.perf_across,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(res.perf_across,1),'FaceColor','b')
plot([0 size(res.perf_across,2)+1],[100 100], '-.', 'color', [.7 .7 .7], 'linewidth',1.5)
plot([0 size(res.perf_across,2)+1],[Pchance Pchance], '-.k', 'linewidth',1.5)
xlabel('Parameter b'),ylabel('Across-session accuracy (%)')
ylim([0 110])
set(gca,'TickDir','out','FontSize',14), box off
title(['History bins: ' num2str(hist_bins) ' - Bin size: ' num2str(td1(1).bin_size)])

subplot(333), hold on
errorbar(mean(res.norm_perf_across,1),std(res.norm_perf_across,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(res.norm_perf_across,1),'FaceColor','g')
plot([0 size(res.norm_perf_across,2)+1],[100 100], 'color', 'k', 'linewidth',.5)
ylabel('Across-session accuracy (%)')
xlabel('Parameter b'),ylim([0 150])
set(gca,'TickDir','out','FontSize',14), box off
title(['History bins: ' num2str(hist_bins) ' - Bin size: ' num2str(td1(1).bin_size)])


subplot(334),hold on
hw = bar(x_hist(1:end-1),hist_within,'histc');
% hn = bar(x_hist(1:end-1),hist_noxval,'histc');
yl = ylim;
set(hw,'FaceColor','k'); alpha(hw,0.7)
%set(hn,'FaceColor','r'); alpha(hn,0.7)
plot([Pchance, Pchance],[yl(1) yl(2)],'-.k', 'linewidth',1.5)
% legend('X-val','Not','Location','NorthWest'), %legend boxoff
xlabel(['Within-session accuracy (%) b = ' num2str(size(perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 x_hist(end-1)])

subplot(335), hold on
ha = bar(x_hist(1:end-1),hist_across,'histc');
yl = ylim;
plot([Pchance, Pchance],[yl(1) yl(2)],'-.k', 'linewidth',1.5)
set(ha,'FaceColor','b');
xlabel(['Across-session accuracy (%) b = ' num2str(size(res.perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 x_hist(end-1)])

subplot(336), hold on
ha = bar(x_hist_norm(1:end-1),hist_norm_across,'histc');
set(ha,'FaceColor','g');
xlabel(['Normalized across accuracy (%) b = ' num2str(size(res.perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
% xlim([0 x_hist_norm(end-1)])
xlim([0 x_hist_norm(end-1)])

subplot(337),hold on
plot(1.25*ones(1,length(res.err_within_xval1_m)),rad2deg(res.err_within_xval1_m),'.','markersize',20,'color',[.7 .7 .7],'linestyle','none')
errorbar(rad2deg(mean(res.err_within_xval1_m,1)),rad2deg(std(res.err_within_xval1_m)),'.k','markersize',30, 'linestyle', 'none', 'linewidth', 1.5)
xlabel('Parameter b'), ylabel('Angular error')
xlim([0 size(perf_within,2)+1])
set(gca,'TickDir','out','FontSize',14), box off
yl = ylim;
ylim([0 floor(yl(2)/10)*10])

subplot(338),hold on
plot(1.25*ones(1,length(res.err_across)),rad2deg(res.err_across),'.','markersize',20,'color','c','linestyle','none')
errorbar(rad2deg(mean(res.err_across,1)),rad2deg(std(res.err_across)),'.b','markersize',30, 'linestyle', 'none', 'linewidth', 1.5)
xlabel('Parameter b'), ylabel('Angular error')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 size(perf_within,2)+1])
yl = ylim;
ylim([0 floor(yl(2)/10)*10])

subplot(339),hold on
plot(x_fit,y_across,'b','linewidth',1.5)
plot(x_fit,y_within2,'color',[.7 .7 .7],'linewidth',1.5)
plot(diff_days,res.perf_across,'.b','markersize',32)
plot(diff_days,res.perf_within_xval2_m,'.','markersize',32,'color',[.7 .7 .7])
ylim([0 100]),xlim([0 max(diff_days)+1])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days since classifier training'), ylabel('Accuracy (%)')
set(gcf,'color','w')
legend('across','within','Location','SouthWest'),legend boxoff
   

ff = '/Users/juangallego/Dropbox/Juan and Matt - Stability latent activity/Results/PMd/';
fn = [td(1).monkey '_' clas_input(1:end-10) '_' method '_' num2str(win_size*1000) 'ms_bin'];

savefig(f1,[ff filesep fn]);
saveas(f1,[ff filesep fn '.png']);

% end



% % --------------------------------------------------------------------------------------------------
% % --------------------------------------------------------------------------------------------------
% % Other plots to understand what's going on
% 
% 
% figure; hold on; cols = parula(9); 
% for t = 1:8
%     for r = 1:15
%         idx = r + (t-1)*15; 
%         plot3( X1(idx,1), X1(idx,2), X1(idx,3), '.', 'color', cols(t,:),'markersize',30); 
%     end 
% end
% 
% % figure; hold on;
% % for t = 1:8
% %     for r = 1:15
% %         idx = r + (t-1)*15; 
% %         plot3(td1(idx).PMd_pca_align(:,1),td1(idx).PMd_pca_align(:,2),td1(idx).PMd_pca_align(:,3), ...
% %             '.', 'color', cols(t,:),'markersize',30); 
% %     end
% % end
% 
% figure; hold on;
% for t = 1:8
%     for r = 1:15
%         idx = r + (t-1)*15; 
%         plot3(td_orig(idx).PMd_pca(:,1),td_orig(idx).PMd_pca(:,2),td_orig(idx).PMd_pca(:,3), ...
%             'color', cols(t,:)); 
%         plot3(td_orig(idx).PMd_pca(end,1),td_orig(idx).PMd_pca(end,2),td_orig(idx).PMd_pca(end,3), ...
%             '.', 'color', cols(t,:),'markersize',30); 
%     end
%     %pause
% end


% clear perf_* err*