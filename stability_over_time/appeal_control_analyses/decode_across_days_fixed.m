%
% Decode from aligned latent signals
%

function results = decode_across_days_fixed( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
in                      = [];
out                     = [];
hist_bins               = 1:5;
n_folds                 = 6; % 'cca' or 'procrustes'
manifold                = [];
mani_dims               = 1:10;
unsort_chs_to_pred      = false;
lag_kin_S1              = 0.050; % in s

if nargin > 1, assignParams(who,params); end % overwrite defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(in), error('Must provide decoder inputs'); end
if isempty(out), error('Must provide decoder output'); end
if isempty(manifold), error('Must provide decoder output'); end
% ToDo: Implement some more checks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STUFF

% get all pairs of sessions
sessions                = unique({td.date});
n_sessions              = length(sessions);

% Sometimes the sessions are not sorted by time --fix that here by
% resorting the trials in master_td
[~, i_sort]         = sort(datenum([sessions]));
if sum( i_sort - 1:length(i_sort) ) > 0
    
    sorted_dates = sort( cell2mat( cellfun(@(x) datenum(x), sessions, 'uni', 0) ) );
    for s = 1:n_sessions
        sessions{s} = datestr(sorted_dates(s),'mm-dd-yyyy');
    end
end

comb_sessions           = nchoosek(1:n_sessions,2);

% Set model parameters
mod_params.model_type   = 'linmodel';
mod_params.out_signals  = out;

% Convert the kinematics to S1 lag to number of bins
lag_kin_S1_bins         = 1:round( lag_kin_S1 / td(1).bin_size );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

disp(['Testing decoder stability using ' in ' inputs']);
disp(['Decoding: ' out]);

c = 0;
for iSess1 = 1:n_sessions
    comb_sessions_temp = comb_sessions(comb_sessions(:,1)==iSess1,:);
    n_comb_sessions         = size(comb_sessions_temp,1);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For all pairs of sessions, test it on session 2
    for iSess2 = 1:n_comb_sessions
        c = c+1;
        
        % get two TD structs, one per session to compare
        [trials2, td2]      = getTDidx(td,'date',sessions{comb_sessions_temp(iSess2,2)});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % train the day 1 decoder
        if iSess1 == 1 && iSess2 == 1
            
            [trials1, td1]      = getTDidx(td,'date',sessions{iSess1});
            
            % a) Set model params inputs
            if max(abs(hist_bins)) ~= 0
                mod_params.in_signals = {[manifold '_align_shift'],1:length(mani_dims)*hist_bins};
            else
                mod_params.in_signals = {[manifold '_align'],1:length(mani_dims)};
            end
            
            % b) "align" dynamics with CCA
            cca_info1 = compDynamics( td, manifold, trials1, trials2, mani_dims );
            
            
            % c) Add latent activity to the tdi structs
            bins_p_trial = size(td1(1).pos,1);
            
            
            for t = 1:length(trials2)
                be      = bins_p_trial*(t-1)+1;
                en      = bins_p_trial*t;
                temp = cca_info1.U(be:en,:);
                td1(t).([manifold '_align' ]) = temp*inv(cca_info1.A) + ...
                        repmat(mean(getSig(td(trials1),manifold),1),size(temp,1),1);
            end
            
            % d) Duplicate and shift to have history
            if max(abs(hist_bins)) ~= 0
                td1_shift         = dupeAndShift(td1,[manifold '_align'],hist_bins,manifold,hist_bins);
            else
                td1_shift = td1;
            end
            
            % past kinematics for S1
            if strcmp( params.manifold(1:end-4), 'S1' ) && ...
                    lag_kin_S1_bins ~= 0 && ...
                    sum( strcmp( params.out, {'pos','vel','acc'} ) ) == 1
                
                % duplicate and shift to get past kinematics
                td1_shift             = dupeAndShift(td1,{out,lag_kin_S1_bins});
                
                % update model outputs to past kinematics
                mod_params.out_signals = {[out '_shift'], [size(td1_shift(1).vel_shift,2)-1 size(td1_shift(1).vel_shift,2)]};
            end
            
            
            % build the Day 1 model
            [td1_shift, mod_info]     = getModel(td1_shift,mod_params);
            
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. PREPARE MODEL INPUTS
        
        switch in
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FOR THE ALIGNED LATENT ACTIVITY
            case 'aligned_data'
                
                
                % a) Set model params inputs
                if max(abs(hist_bins)) ~= 0
                    mod_params.in_signals = {[manifold '_align_shift'],1:length(mani_dims)*hist_bins};
                else
                    mod_params.in_signals = {[manifold '_align'],1:length(mani_dims)};
                end
                
                % b) "align" dynamics with CCA
                cca_info(c) = compDynamics( td, manifold, trials1, trials2, mani_dims );
                
                
                % c) Add latent activity to the tdi structs
                bins_p_trial = size(td1(1).pos,1);
                
                for t = 1:length(trials2)
                    be      = bins_p_trial*(t-1)+1;
                    en      = bins_p_trial*t;
                    % project this aligned data back into space of Day 1
                    temp = cca_info(c).V(be:en,:);
                    td2(t).([manifold '_align' ]) = temp*inv(cca_info(c).A) + ...
                        repmat(mean(getSig(td(trials1),manifold),1),size(temp,1),1);
                end
                
                % d) Duplicate and shift to have history
                if max(abs(hist_bins)) ~= 0
                    td2_shift         = dupeAndShift(td2,[manifold '_align' ],hist_bins,manifold,hist_bins);
                else
                    td2_shift = td2;
                end
                
            otherwise
                error('wrong decoder input');
        end
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1B. IF IT'S FOR S1, WE WANT TO PREDICT PAST KINEMATICS !!!
        
        if strcmp( params.manifold(1:end-4), 'S1' ) && ...
                lag_kin_S1_bins ~= 0 && ...
                sum( strcmp( params.out, {'pos','vel','acc'} ) ) == 1
            
            % duplicate and shift to get past kinematics
            td2_shift             = dupeAndShift(td2,{out,lag_kin_S1_bins});
            
            % update model outputs to past kinematics
            mod_params.out_signals = {[out '_shift'], [size(td1(1).vel_shift,2)-1 size(td1(1).vel_shift,2)]};
        end
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % b) As a reference build a model on day 2. This will serve to
        % normalize the across day predictions obtained by using the decoder
        % trained on day 1 to predict data from day 2
        % Trick: first randomize trial order, because the function that
        % equalizes the number of trials across targets and session sorts them
        % by target
        
        idx_rnd             = randperm(length(td2_shift));
        td2_shift                 = td2_shift(idx_rnd);
        trials_p_fold       = floor(length(td2_shift)/n_folds);
        
        R2xval              = zeros(n_folds,size(td2_shift(1).(out),2));
        
        for f = 1:n_folds
            trials_test     = (f-1)*trials_p_fold+1 : f*trials_p_fold;
            trials_train    = setdiff(1:length(td2_shift),trials_test);
            td_test         = td2_shift(trials_test);
            td_train        = td2_shift(trials_train);
            
            [~, mod_info_xval] = getModel(td_train,mod_params);
            % [R2train(f,:), ~] = testModel(td_train,mod_info_xval);
            
            % compute the cross-validated R2
            [R2xval(f,:), ~] = testModel(td_test,mod_info_xval);
        end
        
        
        % Add XVal results to struct
        res.withinR2_m(c,:) = mean(R2xval,1);
        res.withinR2_sd(c,:) = std(R2xval,0,1);
        res.withinR2_all{c} = R2xval;
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. TEST IT ON DAY 2 -- this is by definition cross-validated
        [R2diff, td2_shift]       = testModel(td2_shift,mod_info);
        
        % Add results to struct
        res.acrossR2(c,:)   = R2diff;
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. CONTROL: compute how well a decoder generalizes without aligning
        % the latent activity
        
        if strcmp(in,'aligned_data')
            
            % a) duplicate and shift --to have history in the model
            if max(abs(hist_bins)) ~= 0
                ctrl_td1            = dupeAndShift(td1,{manifold,hist_bins});
                ctrl_td2            = dupeAndShift(td2,{manifold,hist_bins});
                ctrl_mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
            else
                ctrl_td1 = td1;
                ctrl_td2 = td2;
                ctrl_mod_params.in_signals = {manifold,1:length(mani_dims)};
            end
            
            % b) build model on day 1
            ctrl_mod_params.model_type = mod_params.model_type;
            ctrl_mod_params.out_signals = mod_params.out_signals;
            
            
            [~, ctrl_mod_info]  = getModel(ctrl_td1,ctrl_mod_params);
            
            % c) test it on day 2
            R2ctrl              = testModel(ctrl_td2,ctrl_mod_info);
            res.ctrlR2(c,:)     = R2ctrl;
        end
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5. SEE IF THERE IS A GAIN CHANGE IN THE MODEL
        % -- This could happen if target distance to the center has changed, as
        % it did in Chewie :'-/
        
        Yacross             = get_vars(td2_shift,{out,1:size(td2_shift(1).(out),2)});
        Yhat_across         = get_vars(td2_shift,{'linmodel_default',1:size(td2_shift(1).(out),2)});
        
        res.gain(c,:)       = computeGain(Yacross,Yhat_across);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VARIABLES

% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end
res.diff_days           = diff_days;
res.comb_sessions       = comb_sessions;


% return variable
results                 = res;