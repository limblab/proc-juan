%
%
%

function cc_within = canon_corr_within_task( single_trial_data, dims, varargin )


% -------------------------------------------------------------------------
% read inputs
if nargin >= 3
    shuff           = varargin{1};
else
    shuff           = false;
end
if nargin >= 4
    target          = varargin{2};
else
    target          = 'all_conc'; % by default do all concatenated trials
end
if nargin == 5
    time_window     = varargin{3};
    if size(time_window) ~= [1, 2], error('time_window has wrong dimensions'); end
else
    time_window     = [0 0.5];
end

                        
% -------------------------------------------------------------------------
% get some meta info
nbr_bdfs            = length(single_trial_data);


% -------------------------------------------------------------------------
% prepare data for doing canonical correlation

% 1) equalize trial duration across all tasks
single_trial_data   = equalize_single_trial_dur( single_trial_data, ...
                        'time_win', time_window );
                        

% -------------------------------------------------------------------------
% do

for t = 1:nbr_bdfs
   
    split_data      = split_single_trial_data( single_trial_data{t}, 2, shuff );
    
    % split the data from one task into two
    switch target
        case 'all_conc'
            scores_1 = split_data{1}.target{end}.neural_data.dim_red.scores(:,dims);
            scores_2 = split_data{2}.target{end}.neural_data.dim_red.scores(:,dims);
        otherwise
            scores_1 = split_data{target}.target{end}.neural_data.dim_red.scores(:,dims);
            scores_2 = split_data{target}.target{end}.neural_data.dim_red.scores(:,dims);
    end
            
    % do CCA
   [A, B, r, U, V, stats] = canoncorr( scores_1, scores_2 ); 
   
   
   % store, to return
   cc_within.lin_transform(t).A = A;
   cc_within.lin_transform(t).B = B;
   cc_within.cc(t,:)            = r;
   cc_within.lin_transform(t).U = U;
   cc_within.lin_transform(t).V = V;
   cc_within.stats(t)           = stats;
end