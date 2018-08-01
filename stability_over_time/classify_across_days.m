%
%
%
% ZSCORING DOESN'T WORK WELL


function results = classify_across_days( td, which_type, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
out                 = 'target_direction';
method              = 'Bayes'; % 'NN' 'Bayes'
hist_bins           = 0;
n_folds             = 100; % 'cca' or 'procrustes'
manifold            = [];
mani_dims           = 1:10;
idx_start_align     = {};
idx_end_align       = {};
idx_start_classify  = {};
idx_end_classify    = {};
num_test_trials     = 1; % number per target for within session  CV
if nargin > 1, assignParams(who,params); end % overwrite defaults



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hist_bins < 0, error('Must provide decoder output'); end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

disp(['Testing classifier stability using ' which_type ' inputs']);
disp(['Type: ' method]);
disp(['Classifying: ' out]);


% For all pairs of sessions, build a model on session 1 and test it on
% session 2
for c = 1:n_comb_sessions
    
    disp(['Combination ' num2str(c) '/' num2str(n_comb_sessions)]);
    
    % =====================================================================
    % CLASSIFIER BASED ON ALIGNED LATENT ACTIVITY
    % =====================================================================
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALIGN LATENT ACTIVITY -- HAS TO BE DONE WITH THE ORIGINAL DATA,
    %                   BECAUSE THAT'S WHAT WE USE FOR THE OTHER ANALYSES
    td = td_orig;
    
    switch lower(which_type)
        case {'aligned_data','unaligned_data'}
            % a) get two TD structs, one per session to compare
            trials1         = getTDidx(td,'date',sessions{comb_sessions(c,1)});
            trials2         = getTDidx(td,'date',sessions{comb_sessions(c,2)});
            
            % b) "align" dynamics with CCA
            if strcmpi(which_type,'aligned_data')
                model_in = [manifold '_align'];
                td = trimTD(td,idx_start_align,idx_end_align);
                cca_info(c)     = compDynamics( td, manifold, trials1, trials2, mani_dims );
            else
                model_in = manifold;
            end
            
            % c) Add latent activity to the tdi structs
            td = td_orig;
            td1             = td(trials1);
            td2             = td(trials2);
            
            if strcmpi(which_type,'aligned_data')
                for t = 1:length(td1)
                    temp = td1(t).(manifold);
                    proj1       = temp(:,mani_dims) * cca_info(c).A;
                    td1(t).(model_in) = proj1;
                end
                for t = 1:length(td2)
                    temp = td2(t).(manifold);
                    proj2       = temp(:,mani_dims) * cca_info(c).B;
                    td2(t).(model_in) = proj2;
                end
            end
            
            if hist_bins > 0
                td1      	= dupeAndShift(td1,model_in,hist_bins);
                td2        	= dupeAndShift(td2,model_in,hist_bins);
            end
            
            for trial = 1:length(td1)
                temp = td1(trial).(model_in);
                td1(trial).(model_in) = temp(:,mani_dims);
            end
            for trial = 1:length(td2)
                temp = td2(trial).(model_in);
                td2(trial).(model_in) = temp(:,mani_dims);
            end
            
            td1 = trimTD(td1,idx_start_classify,idx_end_classify);
            td2 = trimTD(td2,idx_start_classify,idx_end_classify);
            
            if hist_bins > 0
                clas_input  = [model_in '_shift'];
            else
                clas_input 	= model_in;
            end
            
            n_clas_vars     = mani_dims;
            
            % define some new parameters to pass throughout the function
            new_params = struct( ...
                'method',method, ...
                'n_folds',n_folds, ...
                'clas_input',clas_input, ...
                'n_clas_vars',n_clas_vars, ...
                'hist_bins',hist_bins, ...
                'num_test_trials',num_test_trials);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CROSS-VALIDATED WITHIN-DAY PREDICTIONS FOR DAY 1
            %       - This is currently done without optimizing the classifier
            
            [~, perf_test1, err_test1] = do_classify_xval(td1,new_params);
            
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TRAIN CLASSIFIER ON DAY 1
            [X1,Y1] = process_td_for_classification(td1,new_params);
            
            % b) Train classifier
            clas = fit_classifier(X1,Y1,new_params);
            
            
            % c) Test performance on training set -- Remember that the optimal
            %       parameters where identified with cross-validation, thus these
            %       predictions are cross-validated
            pred_within = clas.predict(X1)';
            %pred_within     = predict(clas,X1)';
            perf_within     = 100 * ( 1 - sum( angleDiff(pred_within,Y1,true,true) ~= 0 ) / length(pred_within) );
            
            
            % d) Compute classif error (degrees)
            err_within      = mean(angleDiff(Y1,pred_within,true,false));
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TEST ON DAY 2
            
            [X2,Y2] = process_td_for_classification(td2,new_params);
            
            
            % b) Test performance on day 2
            pred_across = clas.predict(X2)';
            
            perf_across     = 100 * ( 1 - sum( angleDiff(pred_across,Y2,true,true) ~= 0 ) / length(pred_across) );
            
            % c) Classifier error
            err_across      = mean(angleDiff(Y2,pred_across,true,false));
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WITHIN-DAY DECODER FOR DAY 2, TO NORMALIZE THE TEST RESULTS
            [~, perf_test2, err_test2] = do_classify_xval(td2,new_params);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Save results
            
            
            % Optimized classifier, not cross-validated on day 1
            res.perf_within_opt_noxval1(c) = perf_within;
            res.err_within_opt_noxval1(c) = err_within;
            
            % Cross-validated on day 1
            res.perf_within_xval1(c,:) = perf_test1;
            res.err_within_xval1(c,:) = err_test1;
            
            % Cross-validated on day 2
            res.perf_within_xval2(c,:) = perf_test2;
            res.err_within_xval2(c,:) = err_test2;
            
            % Across day classifier
            res.perf_across(c) = perf_across;
            res.err_across(c) = err_across;
            
            % For the normalized predictions
            res.norm_perf_across(c) = perf_across/mean(perf_test2)*100;
            
            
            
        case 'spikes'
            td = td_orig;
            
            % a) get two TD structs, one per session to compare
            trials1         = getTDidx(td,'date',sessions{comb_sessions(c,1)});
            trials2         = getTDidx(td,'date',sessions{comb_sessions(c,2)});
            
            td1             = td(trials1);
            td2             = td(trials2);
            
            if hist_bins > 0
                td1      	= dupeAndShift(td1,model_in,hist_bins);
                td2        	= dupeAndShift(td2,model_in,hist_bins);
            end
            
            td1 = trimTD(td1,idx_start_classify,idx_end_classify);
            td2 = trimTD(td2,idx_start_classify,idx_end_classify);
            
            if hist_bins > 0
                clas_input  = [model_in '_shift'];
            else
                clas_input 	= model_in;
            end
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CLASSIFIERS BASED ON SPIKES, WITHIN AND ACROSS DAY
            
            % first get the electrode/unit pairs that are common across days
            sg1 = td1(1).PMd_unit_guide;
            sg2 = td2(1).PMd_unit_guide;
            [~,idx1,idx2] = intersect(sg1,sg2,'rows');
            
            
            %%%% this code removes cells with zero variance in the
            %%%% conditions
            % a bit of a hack to initialize the size of things... just run
            % it once to start
            [X1_spike,Y1_spike] = process_td_for_classification(td1,new_params);
            % check the spikes
            utheta = unique(Y1_spike);
            v = true(length(utheta),size(X1_spike,2));
            
            while any(any(v == true,1))
                [X1_spike,Y1_spike] = process_td_for_classification(td1,new_params);
                [X2_spike,Y2_spike] = process_td_for_classification(td2,new_params);
                
                % check the spikes
                utheta = unique(Y1_spike);
                v1 = zeros(length(utheta),size(X1_spike,2));
                for u = 1:length(utheta)
                    idx = Y1_spike == utheta(u);
                    v1(u,:)  = var(X1_spike(idx,:));
                end
                v2 = zeros(length(utheta),size(X2_spike,2));
                for u = 1:length(utheta)
                    idx = Y2_spike == utheta(u);
                    v2(u,:)  = var(X2_spike(idx,:));
                end
                
                v = v1 == 0 | v2 == 0;
                
                idx = find(any(v,1));
                for u = 1:length(idx)
                    temp = idx(u);
                    if temp > length(idx1), temp = temp - length(idx1); end
                    idx1(temp) = NaN;
                    
                    if temp > length(idx2), temp = temp - length(idx2); end
                    idx2(temp) = NaN;
                end
                
                idx1(isnan(idx1)) = [];
                idx2(isnan(idx2)) = [];
            end
            
            [X1_spike,Y1_spike] = process_td_for_classification(td1,new_params);
            [X2_spike,Y2_spike] = process_td_for_classification(td2,new_params);
            
            % b) Train classifier
            if ~isempty(X1_spike) && ~isempty(X2_spike)
                clas = fit_classifier(X1_spike,Y1_spike,new_params);
                
                % c) Test performance on training set -- Remember that the optimal
                %       parameters where identified with cross-validation, thus these
                %       predictions are cross-validated
                pred_within_spike = clas.predict(X1_spike)';
                perf_within_spike     = 100 * ( 1 - sum( angleDiff(pred_within_spike,Y1_spike,true,true) ~= 0 ) / length(pred_within_spike) );
                
                % d) Compute classif error (degrees)
                err_within_spike      = mean(angleDiff(Y1_spike,pred_within_spike,true,false));
                
                % b) Test performance on day 2
                pred_spike = clas.predict(X2_spike)';
                
                perf_spike     = 100 * ( 1 - sum( angleDiff(pred_spike,Y2_spike,true,true) ~= 0 ) / length(pred_spike) );
                err_spike      = mean(angleDiff(Y2_spike,pred_spike,true,false));
                
                
                % Spike classifier
                res.perf_spike(c) = perf_spike;
                res.err_spike(c) = err_spike;
                % res.norm_perf_spike(c) = perf_spike/mean(perf_test2)*100;
            else
                res.perf_spike(c) = 0;
                res.err_spike(c) = 0;
            end
            
    end
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY STATS AND PLOTS

switch lower(which_type)
    case 'align'
        % Compute mean training performance per day
        res.perf_within_xval1_m  = mean(res.perf_within_xval1,2);
        res.perf_within_xval2_m  = mean(res.perf_within_xval2,2);
        
        res.err_within_xval1_m  = mean(res.err_within_xval1,2);
        res.err_within_xval2_m  = mean(res.err_within_xval2,2);
        
        
        % get number of days between sessions, to summarize the results
        diff_days   = zeros(1,size(comb_sessions,1));
        for c = 1:size(comb_sessions,1)
            diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
        end
        res.diff_days           = diff_days;
        res.comb_sessions       = comb_sessions;
        
    case 'spikes'
        % nothing for now
end

results                 = res;



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y] = process_td_for_classification(td,params)
hist_bins = params.hist_bins;
n_clas_vars = params.n_clas_vars;
clas_input = params.clas_input;

% ii) Prepare the data for the classifier -- Dimensions:
% X: trial x datapoints; Y: trial x target
X = zeros(length(td),length(n_clas_vars)*(hist_bins+1));
for i = 1:length(td)
    
    temp = get_the_signals(td(i),clas_input,n_clas_vars);
    
    X(i,:) = mean(temp,1);
end

Y  = [td.(out)];

X  = zscore(X);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [perf_train, perf_test, err_test] = do_classify_xval(td,params)
n_folds = params.n_folds;
num_test_trials = params.num_test_trials;

% To store the results
perf_train      = zeros(1,n_folds);
perf_test       = zeros(1,n_folds);
err_test       = zeros(1,n_folds);

% b) Do
for f = 1:n_folds
    % a) To cross-validate, shuffle the trials because they are sorted by
    % target location
    [itrain, itest] = get_test_train_trials(td,num_test_trials);
    
    % i) Retrieve the training and testing bins
    
    [Xtrain,Ytrain] = process_td_for_classification(td(itrain),params);
    
    % iii) Train the classifier
    clas_xval = fit_classifier(Xtrain,Ytrain,params);
    
    [Xtest,Ytest] = process_td_for_classification(td(itest),params);
    
    pred_test = clas_xval.predict(Xtest)';
    % Performance (% succesfully classified targets)
    perf_test(f) = 100 * ( 1 - sum( angleDiff(pred_test,Ytest,true,true) ~= 0 ) / length(pred_test) );
    err_test(f) = mean(angleDiff(Ytest,pred_test,true,false));
    
    % v) For reference, compute performance on the training set
    pred_train  = clas_xval.predict(Xtrain)';
    perf_train(f) = 100 * ( 1 - sum( angleDiff(pred_train,Ytrain,true,true) ~= 0 ) / length(pred_train) );
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clas = fit_classifier(X,Y,params)
method = params.method;

switch method
    case 'NN'
        clas = fitcknn(X,Y);
    case 'Bayes'
        clas = fitcnb(X,Y,'ClassNames',unique(Y), ...
            'DistributionNames','normal');
    otherwise
        error('Only NN implemented so far');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = get_the_signals(td,params)
clas_input = params.clas_input;
n_clas_vars = params.n_clas_vars;

if length(td) > 1
    error('currently assumes you pass only one trial.');
end

% see if it's a shift
if ~isempty(regexp(clas_input,'_shift', 'once'))
    % use present bin and shifted
    temp1 = td.(clas_input(1:end-6));
    temp2 = td.(clas_input);
    
    % we need to grab the indices from all the shifted signals
    orig_n = size(temp2,2)/size(temp1,2);
    idx = [];
    for n  = 1:orig_n
        idx = [idx, n_clas_vars + size(temp1,2)*(n-1)];
    end
    
    s = [temp1(:,n_clas_vars), temp2(:,idx) ];
else
    s = td.(clas_input);
    s = s(:,n_clas_vars);
end

end


