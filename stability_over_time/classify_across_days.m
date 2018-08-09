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
idx_start           = {};
idx_end             = {};
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
    % a) get two TD structs, one per session to compare
    trials1         = getTDidx(td,'date',sessions{comb_sessions(c,1)});
    trials2         = getTDidx(td,'date',sessions{comb_sessions(c,2)});
    
    if strcmpi(which_type,'align')
        % b) "align" dynamics with CCA
        td = trimTD(td,idx_start,idx_end);
        cca_info(c)     = compDynamics( td, manifold, trials1, trials2, mani_dims );
    end
    
    % now get the td structs for the two comparison days
    td = td_orig;
    td1             = td(trials1);
    td2             = td(trials2);
    
    switch lower(which_type)
        case 'align'
            % c) Add latent activity to the tdi structs
            for t = 1:length(td1)
                temp = td1(t).(manifold);
                proj1       = temp(:,mani_dims) * cca_info(c).A;
                td1(t).([manifold '_align']) = proj1;
            end
            for t = 1:length(td2)
                temp = td2(t).(manifold);
                proj2       = temp(:,mani_dims) * cca_info(c).B;
                td2(t).([manifold '_align']) = proj2;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Do aligned test and save results
            temp_results = do_manifold_classification_stuff(td1,td2,[manifold '_align'],params);
            
            % Optimized classifier, not cross-validated on day 1
            res.perf_within_opt_noxval1(c) = temp_results.perf_within;
            res.err_within_opt_noxval1(c) = temp_results.err_within;
            
            % Cross-validated on day 1
            res.perf_within_xval1(c,:) = temp_results.perf_test1;
            res.err_within_xval1(c,:) = temp_results.err_test1;
            
            % Cross-validated on day 2
            res.perf_within_xval2(c,:) = temp_results.perf_test2;
            res.err_within_xval2(c,:) = temp_results.err_test2;
            
            % Across day classifier
            res.perf_across(c) = temp_results.perf_across;
            res.err_across(c) = temp_results.err_across;
            
            % the real directions
            res.dir_across(c,:) = temp_results.dir_across';
            res.dir_within(c,:) = temp_results.dir_within';
            res.pred_across(c,:) = temp_results.pred_across;
            res.pred_within(c,:) = temp_results.pred_within;
            
            % For the normalized predictions
            res.norm_perf_across(c) = temp_results.perf_across/mean(temp_results.perf_test2)*100;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Do unaligned control test and save results
            temp_results = do_manifold_classification_stuff(td1,td2,manifold,params);
            
            % Optimized classifier, not cross-validated on day 1
            res.perf_within_opt_noxval1_ctrl(c) = temp_results.perf_within;
            res.err_within_opt_noxval1_ctrl(c) = temp_results.err_within;
            
            % Cross-validated on day 1
            res.perf_within_xval1_ctrl(c,:) = temp_results.perf_test1;
            res.err_within_xval1_ctrl(c,:) = temp_results.err_test1;
            
            % Cross-validated on day 2
            res.perf_within_xval2_ctrl(c,:) = temp_results.perf_test2;
            res.err_within_xval2_ctrl(c,:) = temp_results.err_test2;
            
            % Across day classifier
            res.perf_across_ctrl(c) = temp_results.perf_across;
            res.err_across_ctrl(c) = temp_results.err_across;
            
                        % the real directions
            res.dir_across_ctrl(c,:) = temp_results.dir_across';
            res.dir_within_ctrl(c,:) = temp_results.dir_within';
            res.pred_across_ctrl(c,:) = temp_results.pred_across;
            res.pred_within_ctrl(c,:) = temp_results.pred_within;
            
            % For the normalized predictions
            res.norm_perf_across_ctrl(c) = temp_results.perf_across/mean(temp_results.perf_test2)*100;
            
            
        case 'spikes'
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CLASSIFIERS BASED ON SPIKES
            temp_results = do_manifold_classification_stuff(td1,td2,[manifold(1:end-4) '_spikes'],params);
            
            % Optimized classifier, not cross-validated on day 1
            res.perf_within_opt_noxval1_spike(c) = temp_results.perf_within;
            res.err_within_opt_noxval1_spike(c) = temp_results.err_within;
            
            % Cross-validated on day 1
            res.perf_within_xval1_spike(c,:) = temp_results.perf_test1;
            res.err_within_xval1_spike(c,:) = temp_results.err_test1;
            
            % Cross-validated on day 2
            res.perf_within_xval2_spike(c,:) = temp_results.perf_test2;
            res.err_within_xval2_spike(c,:) = temp_results.err_test2;
            
            % Across day classifier
            res.perf_across_spike(c) = temp_results.perf_across;
            res.err_across_spike(c) = temp_results.err_across;
            
            % the real directions
            res.dir_across_spike(c,:) = temp_results.dir_across';
            res.dir_within_spike(c,:) = temp_results.dir_within';
            res.pred_across_spike(c,:) = temp_results.pred_across;
            res.pred_within_spike(c,:) = temp_results.pred_within;
            
            % For the normalized predictions
            res.norm_perf_across_spike(c) = temp_results.perf_across/mean(temp_results.perf_test2)*100;
            
    end
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY STATS AND PLOTS
if strcmpi(which_type,'align')
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
end
results                 = res;



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y] = process_td_for_classification(td,idx,params)
out = params.out;
hist_bins = params.hist_bins;
params.idx = idx;

% ii) Prepare the data for the classifier -- Dimensions:
% X: trial x datapoints; Y: trial x target
X = zeros(length(td),length(idx)*(hist_bins+1));
for i = 1:length(td)
    
    temp = get_the_signals(td(i),params);
    
    X(i,:) = mean(temp,1);
end

Y  = [td.(out)];

X  = zscore(X);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [perf_train, perf_test, err_test] = do_classify_xval(td,idx,params)
n_folds = params.n_folds;
num_test_trials = params.num_test_trials;
clas_input = params.clas_input;

params.idx = idx;

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
    if ~isempty(regexp(clas_input,'_spikes', 'once'))
        temp_idx = remove_zero_var_neurons(td(itrain),params);
    else
        temp_idx = idx;
    end
    params.idx = temp_idx;
    [Xtrain,Ytrain] = process_td_for_classification(td(itrain),temp_idx,params);
    
    % iii) Train the classifier
    if ~isempty(Xtrain)
        clas_xval = fit_classifier(Xtrain,Ytrain,params);
        
        [Xtest,Ytest] = process_td_for_classification(td(itest),temp_idx,params);
        
        pred_test = clas_xval.predict(Xtest)';
        % Performance (% succesfully classified targets)
        perf_test(f) = 100 * ( 1 - sum( angleDiff(pred_test,Ytest,true,true) ~= 0 ) / length(pred_test) );
        err_test(f) = mean(angleDiff(Ytest,pred_test,true,false));
        
        % v) For reference, compute performance on the training set
        pred_train  = clas_xval.predict(Xtrain)';
        perf_train(f) = 100 * ( 1 - sum( angleDiff(pred_train,Ytrain,true,true) ~= 0 ) / length(pred_train) );
    else
        perf_train = NaN;
        perf_test = NaN;
        err_test = NaN;
    end
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
idx = params.idx;

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
        idx = [idx, idx + size(temp1,2)*(n-1)];
    end
    
    s = [temp1(:,idx), temp2(:,idx) ];
else
    s = td.(clas_input);
    s = s(:,idx);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = do_manifold_classification_stuff(td1,td2,model_in,params)
hist_bins = params.hist_bins;
mani_dims = params.mani_dims;
idx_start_classify = params.idx_start_classify;
idx_end_classify = params.idx_end_classify;
out = params.out;
n_folds = params.n_folds;
method = params.method;
num_test_trials = params.num_test_trials;


if hist_bins > 0
    td1      	= dupeAndShift(td1,model_in,hist_bins);
    td2        	= dupeAndShift(td2,model_in,hist_bins);
end

if ~isempty(regexp(model_in,'_spikes','ONCE'))
    % do nothing
else
    % take the manifold dimensions requested
    for trial = 1:length(td1)
        temp = td1(trial).(model_in);
        td1(trial).(model_in) = temp(:,mani_dims);
    end
    for trial = 1:length(td2)
        temp = td2(trial).(model_in);
        td2(trial).(model_in) = temp(:,mani_dims);
    end
end

td1 = trimTD(td1,idx_start_classify,idx_end_classify);
td2 = trimTD(td2,idx_start_classify,idx_end_classify);

if hist_bins > 0
    clas_input  = [model_in '_shift'];
else
    clas_input 	= model_in;
end

% define some new parameters to pass throughout the function
new_params = struct( ...
    'out',out, ...
    'method',method, ...
    'n_folds',n_folds, ...
    'clas_input',clas_input, ...
    'hist_bins',hist_bins, ...
    'num_test_trials',num_test_trials);

if ~isempty(regexp(model_in,'_spikes','ONCE'))
    [idx1, idx2] = remove_zero_var_neurons(td1,td2,new_params);
else
    idx1     = mani_dims;
    idx2     = mani_dims;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROSS-VALIDATED WITHIN-DAY PREDICTIONS FOR DAY 1
%       - This is currently done without optimizing the classifier
if isempty(regexp(model_in,'_spikes','ONCE'))
    [~, res.perf_test1, res.err_test1] = do_classify_xval(td1,mani_dims,new_params);
else
    res.perf_test1 = NaN;
    res.err_test1 = NaN;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAIN CLASSIFIER ON DAY 1 AND DAY 2
[X1,Y1] = process_td_for_classification(td1,idx1,new_params);
[X2,Y2] = process_td_for_classification(td2,idx2,new_params);

% b) Train classifier
if ~isempty(X1) && ~isempty(X2)
    % This is a hacky way of getting some CV performance for the confusion
    % matrix
% % %     % first do some cross validation
% % %     [temp_pred, temp_dir] = deal([]);
% % %     for i = 1:23
% % %         [fuckme, fuckyou] = get_test_train_trials(td2,1);
% % %         clas = fit_classifier(X2(fuckme,:),Y2(fuckme),new_params);
% % %         temp_pred = cat(2,temp_pred,clas.predict(X2(fuckyou,:))');
% % %         temp_dir =  cat(2,temp_dir,Y2(fuckyou));
% % %     end
% % %     res.pred_within = temp_pred;
% % %     res.dir_within = temp_dir;
    
    % now train it on day 1 normally
    clas = fit_classifier(X1,Y1,new_params);
    
    % c) Test performance on training set
    res.pred_within = clas.predict(X1)';
    res.dir_within = Y1;
    
    %pred_within     = predict(clas,X1)';
    res.perf_within     = 100 * ( 1 - sum( angleDiff(res.pred_within,res.dir_within,true,true) ~= 0 ) / length(res.pred_within) );
    
    % d) Compute classif error (degrees)
    res.err_within      = mean(angleDiff(res.dir_within,res.pred_within,true,false));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEST ON DAY 2
    
    % b) Test performance on day 2
    res.pred_across = clas.predict(X2)';
    res.perf_across     = 100 * ( 1 - sum( angleDiff(res.pred_across,Y2,true,true) ~= 0 ) / length(res.pred_across) );
    res.dir_across = Y2;
    
    % c) Classifier error
    res.err_across      = mean(angleDiff(Y2,res.pred_across,true,false));
else
    res.pred_within = NaN;
    res.perf_within = NaN;
    res.err_within = NaN;
    res.pred_across = NaN;
    res.perf_across = NaN;
    res.err_across = NaN;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WITHIN-DAY DECODER FOR DAY 2, TO NORMALIZE THE TEST RESULTS
[~, res.perf_test2, res.err_test2] = do_classify_xval(td2,idx2,new_params);


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx1, idx2] = remove_zero_var_neurons(varargin)

if nargin == 3
    td1 = varargin{1};
    td2 = varargin{2};
    new_params = varargin{3};
    
    % first get the electrode/unit pairs that are common across days
    sg1 = td1(1).PMd_unit_guide;
    sg2 = td2(1).PMd_unit_guide;
    [~,idx1,idx2] = intersect(sg1,sg2,'rows');
    
    %%%% this code removes cells with zero variance in the
    %%%% conditions
    % a bit of a hack to initialize the size of things... just run
    % it once to start
    [X1,Y1] = process_td_for_classification(td1,idx1,new_params);
    
    % check the spikes
    utheta = unique(Y1);
    v = true(length(utheta),size(X1,2));
    
    while any(any(v == true,1))
        [X1,Y1] = process_td_for_classification(td1,idx1,new_params);
        [X2,Y2] = process_td_for_classification(td2,idx2,new_params);
        
        % check the spikes
        utheta = unique(Y1);
        v1 = zeros(length(utheta),size(X1,2));
        for u = 1:length(utheta)
            v1(u,:)  = var(X1(Y1 == utheta(u),:));
        end
        v2 = zeros(length(utheta),size(X2,2));
        for u = 1:length(utheta)
            v2(u,:)  = var(X2(Y2 == utheta(u),:));
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
    
elseif nargin == 2
    td1 = varargin{1};
    new_params = varargin{2};
    idx2 = [];
    
    % first get the electrode/unit pairs that are common across days
    idx1 = (1:size(td1(1).PMd_unit_guide,1))';
    
    
    %%%% this code removes cells with zero variance in the
    %%%% conditions
    % a bit of a hack to initialize the size of things... just run
    % it once to start
    [X1,Y1] = process_td_for_classification(td1,idx1,new_params);
    
    % check the spikes
    utheta = unique(Y1);
    v = true(length(utheta),size(X1,2));
    
    while any(any(v == true,1))
        [X1,Y1] = process_td_for_classification(td1,idx1,new_params);
        
        % check the spikes
        utheta = unique(Y1);
        v1 = zeros(length(utheta),size(X1,2));
        for u = 1:length(utheta)
            v1(u,:)  = var(X1(Y1 == utheta(u),:));
        end
        
        v = v1 == 0;
        
        idx = find(any(v,1));
        for u = 1:length(idx)
            temp = idx(u);
            if temp > length(idx1), temp = temp - length(idx1); end
            idx1(temp) = NaN;
        end
        
        idx1(isnan(idx1)) = [];
    end
    
end


end



