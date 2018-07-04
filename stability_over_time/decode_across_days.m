%
% Decode from aligned latent signals
%

function results = decode_across_days( td, params ) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
in                      = [];
out                     = [];
hist_bins               = 5;
n_folds                 = 6; % 'cca' or 'procrustes'
manifold                = [];
mani_dims               = 1:10;


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
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);

% Set model parameters
mod_params.model_type   = 'linmodel';
mod_params.out_signals  = out;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

disp('Testing decoder stability');
disp(['Decoding: ' out]);


% For all pairs of sessions, build a model on session 1 and test it on
% session 2
for c = 1:n_comb_sessions

    
    % get two TD structs, one per session to compare
    [trials1, td1]      = getTDidx(td,'date',sessions{comb_sessions(c,1)});
    [trials2, td2]      = getTDidx(td,'date',sessions{comb_sessions(c,2)});


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. PREPARE MODEL INPUTS
    
    switch in
    
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOR THE ALIGNED LATENT ACTIVITY
        case 'aligned_data'
            
            
            % a) Set model params inputs
            mod_params.in_signals = {['align_' manifold '_shift'],1:length(mani_dims)*hist_bins};


            % b) "align" dynamics with CCA
            cca_info(c) = compDynamics( td, manifold, trials1, trials2, mani_dims );


            % c) Add latent activity to the tdi structs 
            bins_p_trial = size(td1(1).pos,1);

            for t = 1:length(trials1)
                be      = bins_p_trial*(t-1)+1;
                en      = bins_p_trial*t;        
                td1(t).(['align_' manifold]) = cca_info(c).U(be:en,:);
                td2(t).(['align_' manifold]) = cca_info(c).V(be:en,:);
            end        
            
            % d) Duplicate and shift to have history
            td1         = dupeAndShift(td1,{['align_' manifold],hist_bins});
            td2         = dupeAndShift(td2,{['align_' manifold],hist_bins});
            
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOR THE UNALIGNED LATENT ACTIVITY  
        case 'raw_data'
            
            % a) Set model params inputs
            mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
            
            % b) Duplicate and shift to have history
            td1         = dupeAndShift(td1,{manifold,hist_bins});
            td2         = dupeAndShift(td2,{manifold,hist_bins});
        
        otherwise
            error('wrong decoder input');
    end
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. BUILD A MODEL ON DAY 1
    
    % a) Build it!
    [td1, mod_info]     = getModel(td1,mod_params);

    % b) Cross-validate the predictions
    % Trick: first randomize trial order, because the function that
    % equalizes the number of trials across targets and session sorts them
    % by target
    
    td1 = td1(randperm(length(td1)));
    trials_p_fold = floor(length(td1)/n_folds);

    R2xval              = zeros(n_folds,size(td1(1).(out),2));
    
    for f = 1:n_folds
        trials_test     = (f-1)*trials_p_fold+1 : f*trials_p_fold;
        trials_train    = setdiff(1:length(td1),trials_test);
        td_test         = td1(trials_test);
        td_train        = td1(trials_train);
        
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
    
    % Do it!
    [R2diff, td2]       = testModel(td2,mod_info);
    
    % Add results to struct
    res.acrossR2(c,:)   = R2diff;

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. CONTROL: compute how well a decoder generalizes without aligning
    
    % HERE HERE HERE
    
    % duplicate and shift --to have history in the model
    ctrl_td1            = dupeAndShift(td1,{manifold,hist_bins});
    ctrl_td2            = dupeAndShift(td2,{manifold,hist_bins});

    % build model on day 1
    ctrl_mod_params.model_type = mod_params.model_type;
    ctrl_mod_params.out_signals = mod_params.out_signals;
    ctrl_mod_params.in_signals = {[manifold '_shift'],1:length(mani_dims)*hist_bins};
    
    [~, ctrl_mod_info]  = getModel(ctrl_td1,ctrl_mod_params);

    % test it on day 2
    R2ctrl              = testModel(ctrl_td2,ctrl_mod_info);
    res.ctrlR2(c,:)     = R2ctrl;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. SEE IF THERE IS A GAIN CHANGE IN THE MODEL
    % -- This could happen if target distance to the center has changed, as
    % it did in Chewie :'-/
    
    Yacross             = get_vars(td2,{out,1:size(td2(1).(out),2)});
    Yhat_across         = get_vars(td2,{'linmodel_default',1:size(td2(1).(out),2)});
    
    res.gain(c,:)       = computeGain(Yacross,Yhat_across);
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