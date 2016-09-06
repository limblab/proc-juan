


function can_corrs = canon_corr_all_manifolds( single_trial_data, ...
                            dims, labels, varargin )

                        
% -------------------------------------------------------------------------
% read inputs
if nargin >= 4
    target          = varargin{1};
else
    target          = 'all_conc'; % by default do all concatenated trials
end
if nargin == 5
    time_window     = varargin{2};
    if size(time_window) ~= [1, 2], error('time_window has wrong dimensions'); end
else
    time_window     = [0 0.5];
end


                        
% -------------------------------------------------------------------------
% get some meta info
nbr_bdfs            = length(single_trial_data);

comb_bdfs           = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs       = size(comb_bdfs,1);



% -------------------------------------------------------------------------
% prepare data for doing canonical correlation

% 1) equalize trial duration across all tasks
single_trial_data   = equalize_single_trial_dur( single_trial_data, ...
                        'time_win', time_window );

% 2) equalize number of trials for all targets of a given task
for i = 1:nbr_bdfs
    single_trial_data{i} = equalize_nbr_trials_p_target( single_trial_data{i} );
end

% 3) equalize number of trials across tasks
single_trial_data   = equalize_nbr_trials_across_tasks( single_trial_data, target );
    

% -------------------------------------------------------------------------
% do

% ToDo: init can_corrs matrix

for p = 1:nbr_comb_bdfs
    
    % get projected, aligned in time scores for this task
    switch target
        case 'all_conc'
            scores_1 = single_trial_data{comb_bdfs(p,1)}.target{end}.neural_data.dim_red.scores(:,dims);
            scores_2 = single_trial_data{comb_bdfs(p,2)}.target{end}.neural_data.dim_red.scores(:,dims);
        otherwise
            scores_1 = single_trial_data{comb_bdfs(p,1)}.target{target}.neural_data.dim_red.scores(:,dims);
            scores_2 = single_trial_data{comb_bdfs(p,2)}.target{target}.neural_data.dim_red.scores(:,dims);            
    end
    
    % do CCA
   [A, B, r, U, V, stats] = canoncorr( scores_1, scores_2 ); 
   
   
   % store, to return
   can_corrs.lin_transform(p).A = A;
   can_corrs.lin_transform(p).B = B;
   can_corrs.cc(p,:)            = r;
   can_corrs.lin_transform(p).U = U;
   can_corrs.lin_transform(p).V = V;
   can_corrs.stats(p)           = stats;
end