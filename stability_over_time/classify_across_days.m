%
%
%
% ZSCORING DOESN'T WORK WELL


function results = classify_across_days( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
in                      = [];
out                     = [];
method                  = 'Bayes'; % 'NN' 'Bayes'
% win                     = [];
hist_bins               = [];
win_size                = [];
n_folds                 = 5; % 'cca' or 'procrustes'
manifold                = [];
mani_dims               = 1:6;
zsc                     = false;


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
sessions                = unique({td.date});
n_sessions              = length(sessions);
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

disp(['Testing classifier stability using ' in ' inputs']);
disp(['Type: ' method]);
disp(['Classifying: ' out]);


% PLAYING WITH THINGS -- NEED TO BE MOVED TO PARAMETERS WHEN THIS IS A REAL
% FUNCTION
td_orig                 = td;

% downsample again
%win_size                = 0.150; %0.150;
down_rate               = round(win_size/td_orig(1).bin_size);

% Only keep the last 'win_size' s of the trial
td                      = trimTD(td_orig,{'end',-down_rate},{'end',0});

% Make one big bin (a la Santhanam et al JNP 2007)
td                      = binTD(td,down_rate);


% % DO FOR ALL POSSIBLE HISTORY BINS --NEEDS TO BE DELETED WHEN FINALIZED
% for h = 1:max(size(td(1).pos,1)-2,1)
%     
%     hist_bins = h;

% For all pairs of sessions, build a model on session 1 and test it on
% session 2
for c = 1:n_comb_sessions


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLASSIFIER BASED ON LATENT ACTIVITY

    % a) get two TD structs, one per session to compare
    [trials1, td1]      = getTDidx(td,'date',sessions{comb_sessions(c,1)});
    [trials2, td2]      = getTDidx(td,'date',sessions{comb_sessions(c,2)});


    % b) "align" dynamics with CCA
    cca_info(c)         = compDynamics( td, manifold, trials1, trials2, mani_dims );


    % c) Add latent activity to the tdi structs
    bins_p_trial        = size(td1(1).pos,1);

    for t = 1:length(trials1)
        be              = bins_p_trial*(t-1)+1;
        en              = bins_p_trial*t;
        td1(t).([manifold '_align']) = cca_info(c).U(be:en,:);
        td2(t).([manifold '_align']) = cca_info(c).V(be:en,:);
    end

    % d) Duplicate and shift to have history
    if hist_bins > 0
        td1             = dupeAndShift(td1,{[manifold '_align'],hist_bins});
        td2             = dupeAndShift(td2,{[manifold '_align'],hist_bins});
    end


    % ---------------------------------------------------------------------
    % CLASSIFY ON DAY 1

    if hist_bins > 0
        clas_input      = [manifold '_align_shift'];
    else
        clas_input      = [manifold '_align'];
    end

    % a) To cross-validate, shuffle the trials because they are sorted by
    % target location
    ishuf               = randperm(length(td1));
    bins_p_fold         = floor(length(td1)/n_folds);


    % TEMP: Do for each bin (after dupe and shifting)
    n_clas_vars         = size(td1(1).(clas_input),2);


    % LOOP OVER BINS IN THE WINDOW
    for b = 1:size(td1(1).(clas_input),1)

        % -----------------------------------------------------------------
        % Cross-validated WITHIN-DAY performance

        perf_train      = zeros(1,n_folds);
        perf_test       = zeros(1,n_folds);        

        for f = 1:n_folds

            % a) Retrieve training and testing bins
            itest       = ishuf( (f-1)*bins_p_fold + 1 : f*bins_p_fold );
            itrain      = ishuf(~ismember(ishuf,itest));

            % b) Prepare data for classifier -- X has to be a trial x
            % datapoints matrix; Y a trial x 1 matrix
            Xtrain      = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td1(itrain)', 'uniformoutput', false ) );
            Ytrain      = [td1(itrain).target_direction];

            % c) Train classifier
            if zsc
                Xtrain  = zscore(Xtrain);
            end
            
            switch method
                case 'NN'
                    clas_within = fitcknn(Xtrain,Ytrain);
                case 'Bayes'
%                     clas_within = fitcnb(Xtrain,Ytrain,'DistributionNames','kernel', ...
%                                                         'Kernel','triangle', ...
%                                                         'Support','unbounded', ...
%                                                         'OptimizeHyperparameters',{'Width'});
%                     close all;
                    clas_within = fitcnb(Xtrain,Ytrain);
                otherwise
                    error('Only NN implemented so far');
            end

            % Test performance on training set -- While developing the
            % classifier -- DELETE DELETE -- DELETE DELETE 
            pred_train  = predict(clas_within,Xtrain)';
            perf_train(f) = ( length(pred_train) - sum( ( pred_train - Ytrain ~= 0 ) ) ) / length(pred_train)*100;

            % d) Test on testing set (derrr) on the same day
            Xtest       = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td1(itest)', 'uniformoutput', false ) );            
            Ytest       = [td1(itest).target_direction];

            if zsc
                Xtest   = zscore(Xtest);
%                Ytest   = zscore(Ytest);
            end
            
            pred_test   = predict(clas_within,Xtest)';

            % Performance (% succesfully classified targets)
            perf_test(f) = ( length(pred_test) - sum( ( pred_test - Ytest ~= 0 ) ) ) / length(pred_test)*100;

            % Classif error (degrees)
            err_test(f) = mean(abs(Ytest-pred_test));
        end
        perf_within_m(c,b) = mean(perf_test);
        err_within_m(c,b) = mean(err_test);


        % -----------------------------------------------------------------
        % ACROSS-DAY PERFORMANCE

        % a) Build within day classifier with all the data
        X1              = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                                td1', 'uniformoutput', false ) );            
        Y1              = [td1.target_direction];

        if zsc
            X1          = zscore(X1);
        end
        
        switch method
            case 'NN'
                clas    = fitcknn(X1,Y1);
            case 'Bayes'
%                 clas    = fitcnb(X1,Y1,'DistributionNames','kernel', ...
%                                         'Kernel','triangle', ...
%                                         'Support','unbounded', ...
%                                         'OptimizeHyperparameters',{'Width'});
%                 close all;
                clas    = fitcnb(X1,Y1);
            otherwise
                error('Only NN implemented so far');
        end

        % see how well it predicts the whole set -- DELETE DELETE 
        pred1           = predict(clas,X1)';
        perf1(c)        = ( length(pred1) - sum( ( pred1 - Y1 ~= 0 ) ) ) / length(pred1)*100;

        % b) Test on day 2
        X2              = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                                td2', 'uniformoutput', false ) );            
        Y2              = [td2.target_direction];

        if zsc
            X2          = zscore(X2);
        end
        
        pred2           = predict(clas,X2)';
        
        
        % c) Train within-session decoder on day 2
        perf_test2      = zeros(1,n_folds);
        
        for f = 1:n_folds
            
            % a) Retrieve training and testing bins
            itest       = ishuf( (f-1)*bins_p_fold + 1 : f*bins_p_fold );
            itrain      = ishuf(~ismember(ishuf,itest));

            % b) Prepare data for classifier -- X has to be a trial x
            % datapoints matrix; Y a trial x 1 matrix
            Xtrain2     = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td1(itrain)', 'uniformoutput', false ) );
            Ytrain2     = [td1(itrain).target_direction];

            % c) Train classifier
            if zsc
                Xtrain2 = zscore(Xtrain2);
            end
            
            switch method
                case 'NN'
                    clas_within = fitcknn(Xtrain2,Ytrain2);
                case 'Bayes'
%                     clas_within = fitcnb(Xtrain,Ytrain,'DistributionNames','kernel', ...
%                                                         'Kernel','triangle', ...
%                                                         'Support','unbounded', ...
%                                                         'OptimizeHyperparameters',{'Width'});
%                     close all;
                    clas_within = fitcnb(Xtrain2,Ytrain2);
                otherwise
                    error('Only NN implemented so far');
            end
            
            % d) Test on testing set (derrr) on the same day
            Xtest2      = cell2mat( arrayfun( @(x) x.(clas_input)(b,1:n_clas_vars), ...
                            td1(itest)', 'uniformoutput', false ) );            
            Ytest2      = [td1(itest).target_direction];
            
            if zsc
                Xtest2  = zscore(Xtest2);
%                Ytest   = zscore(Ytest);
            end
            
            pred_test2  = predict(clas_within,Xtest2)';

            % Performance (% succesfully classified targets)
            perf_test2(f) = ( length(pred_test2) - sum( ( pred_test2 - Ytest2 ~= 0 ) ) ) / length(pred_test2)*100;
        end
        
        % TO DELETE TO DELETE TO DELETE TO DELETE TO DELETE TO DELETE TO DELETE TO DELETE
        perf_within_m2(c,b) = mean(perf_test2);
        

        % Performance (% succesfully classified targets)          
        perf_across(c,b) = ( length(pred2) - sum( ( pred2 - Y2 ~= 0 ) ) ) / length(pred2)*100;

        % Classif error (degrees)
        err_across_m(c,b) = mean(abs(Y2-pred2));
        
        % For the normalized predictions
        norm_perf_across(c,b) = perf_across(c,b)/perf_within_m2(c)*100;
    end
end



% Distirbution predictions 
hist_res                = 5;
x_hist                  = 0:hist_res:100+hist_res;
hist_within             = histcounts(perf_within_m(:,end),x_hist)/numel(perf_within_m)*100;
hist_across             = histcounts(perf_across(:,end),x_hist)/numel(perf_across)*100;

x_hist_norm             = 0:hist_res:100+hist_res+50;
hist_norm_across        = histcounts(norm_perf_across(:,end),x_hist_norm)/numel(norm_perf_across)*100;

% not xvalidated within session
hist_noxval             = histcounts(perf1,x_hist)/length(perf1)*100;

Pchance                 = 100/numel(unique([td1.target_direction]));



figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])    
subplot(331), hold on
errorbar(mean(perf_within_m,1),std(perf_within_m,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(perf_within_m,1),'FaceColor',[.7 .7 .7])
plot([0 size(perf_within_m,2)+1],[100 100], '-.', 'color', [.7 .7 .7], 'linewidth',1.5)
plot([0 size(perf_within_m,2)+1],[Pchance Pchance], '-.k', 'linewidth',1.5)
ylim([0 110])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Parameter b'),ylabel('Within-session accuracy (%)')
title(['Classifier: ' method])

subplot(332), hold on
errorbar(mean(perf_across,1),std(perf_across,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(perf_across,1),'FaceColor','b')
plot([0 size(perf_within_m,2)+1],[100 100], '-.', 'color', [.7 .7 .7], 'linewidth',1.5)
plot([0 size(perf_within_m,2)+1],[Pchance Pchance], '-.k', 'linewidth',1.5)
xlabel('Parameter b'),ylabel('Across-session accuracy (%)')
ylim([0 110])
set(gca,'TickDir','out','FontSize',14), box off
title(['History bins: ' num2str(hist_bins) ' - Bin size: ' num2str(td1(1).bin_size)])

subplot(333), hold on
errorbar(mean(norm_perf_across,1),std(norm_perf_across,1),'k','marker','none', 'linestyle', 'none', 'linewidth', 1.5)
bar(mean(norm_perf_across,1),'FaceColor','g')
plot([0 size(norm_perf_across,2)+1],[1 1], 'color', 'k', 'linewidth',.5)
ylabel('Across-session accuracy (%)')
xlabel('Parameter b'),ylim([0 150])
set(gca,'TickDir','out','FontSize',14), box off
title(['History bins: ' num2str(hist_bins) ' - Bin size: ' num2str(td1(1).bin_size)])


subplot(334),hold on
hw = bar(x_hist(1:end-1),hist_within,'histc');
hn = bar(x_hist(1:end-1),hist_noxval,'histc');
yl = ylim;
set(hw,'FaceColor','k'); alpha(hw,0.7)
set(hn,'FaceColor','r'); alpha(hn,0.7)
plot([Pchance, Pchance],[yl(1) yl(2)],'-.k', 'linewidth',1.5)
legend('X-val','Not','Location','NorthWest'), %legend boxoff
xlabel(['Within-session accuracy (%) b = ' num2str(size(perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 x_hist(end-1)])

subplot(335), hold on
ha = bar(x_hist(1:end-1),hist_across,'histc');
yl = ylim;
plot([Pchance, Pchance],[yl(1) yl(2)],'-.k', 'linewidth',1.5)
set(ha,'FaceColor','b');
xlabel(['Across-session accuracy (%) b = ' num2str(size(perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 x_hist(end-1)])

subplot(336), hold on
ha = bar(x_hist_norm(1:end-1),hist_norm_across,'histc');
set(ha,'FaceColor','g');
xlabel(['Normalized across accuracy (%) b = ' num2str(size(perf_across,2))]), ylabel('Sessions (%)')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 x_hist_norm(end-1)])


subplot(337),hold on
plot(1.25*ones(1,length(err_within_m)),rad2deg(err_within_m),'.','markersize',20,'color',[.7 .7 .7],'linestyle','none')
errorbar(rad2deg(mean(err_within_m,1)),rad2deg(std(err_within_m)),'.k','markersize',30, 'linestyle', 'none', 'linewidth', 1.5)
xlabel('Parameter b'), ylabel('Angular error')
xlim([0 size(perf_within_m,2)+1])
set(gca,'TickDir','out','FontSize',14), box off
yl = ylim;
ylim([0 floor(yl(2)/10)*10])

subplot(338),hold on
plot(1.25*ones(1,length(err_across_m)),rad2deg(err_within_m),'.','markersize',20,'color','c','linestyle','none')
errorbar(rad2deg(mean(err_across_m,1)),rad2deg(std(err_across_m)),'.b','markersize',30, 'linestyle', 'none', 'linewidth', 1.5)
xlabel('Parameter b'), ylabel('Angular error')
set(gca,'TickDir','out','FontSize',14), box off
xlim([0 size(perf_within_m,2)+1])
yl = ylim;
ylim([0 floor(yl(2)/10)*10])

set(gcf,'color','w')

    
    
% end



% RETURN VARS
results.perf_within     = perf_within_m;
results.perf_across     = perf_across;



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
% Other plots to understand what's going on


figure; hold on; cols = parula(9); 
for t = 1:8
    for r = 1:15
        idx = r + (t-1)*15; 
        plot3( X1(idx,1), X1(idx,2), X1(idx,3), '.', 'color', cols(t,:),'markersize',30); 
    end 
end

% figure; hold on;
% for t = 1:8
%     for r = 1:15
%         idx = r + (t-1)*15; 
%         plot3(td1(idx).PMd_pca_align(:,1),td1(idx).PMd_pca_align(:,2),td1(idx).PMd_pca_align(:,3), ...
%             '.', 'color', cols(t,:),'markersize',30); 
%     end
% end

figure; hold on;
for t = 1:8
    for r = 1:15
        idx = r + (t-1)*15; 
        plot3(td_orig(idx).PMd_pca(:,1),td_orig(idx).PMd_pca(:,2),td_orig(idx).PMd_pca(:,3), ...
            'color', cols(t,:)); 
        plot3(td_orig(idx).PMd_pca(end,1),td_orig(idx).PMd_pca(end,2),td_orig(idx).PMd_pca(end,3), ...
            '.', 'color', cols(t,:),'markersize',30); 
    end
    %pause
end


clear perf_* err*