%
% Generate confidence level for CCA, by shuffling the neural units' weights
% onto the neural modes. 
%
%   function signif_boots = bootstrap_weights_canon_corr( single_trial_data, ...
%                           dims, varargin )
%
% Inputs (opt)      : [default]
%   single_trial_data : single trial data for all the tasks (cell array)
%   dim             : manifold dimension
%   (target)        : ['all_conc'] targets that will be compared
%   (time_win)      : [0 0.5] window to calculate CC (s)
%   (nbr_shuffles)  : [1,000] nbr of times it will be repeated
%   (perc_signif)   : [99] percentile that will be given as significance
%                       threshold
%
% Output:
%   signif_boots    : 1-by-n matrix with the CC for each variable that is
%                       obtained by chance
%
%

function signif_boots = bootstrap_weights_canon_corr( single_trial_data, ...
                            dims, varargin )


% -------------------------------------------------------------------------
% parameters

% read inputs
if nargin >= 3
    target          = varargin{1};
else
    target          = 'all_conc'; % by default do all concatenated trials
end

if nargin >= 4
    time_window     = varargin{3};
    if size(time_window) ~= [1, 2], error('time_window has wrong dimensions'); end
else
    time_window     = [0 0.5];
end

if nargin >= 5
    nbr_shuffles    = varargin{2};
end
if nargin == 6
    perc_signif     = varargin{3};
end

% if not passed, use defaults 
if ~exist('nbr_shuffles','var')
    nbr_shuffles    = 1000; 
end
if ~exist('perc_signif','var')
    perc_signif     = 99;
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



% Do for each pair of tasks
for p = 1:nbr_comb_bdfs
    
    cc                  = zeros(length(dims),nbr_shuffles);
    r                   = zeros(length(dims),nbr_shuffles);

    
    % Do for nbr_shuffles times
    for i = 1:nbr_shuffles

        % get projected, aligned in time scores for task #1 and smoothed FRs
        % for task 2, that will then be projected onto shuffled modes
        switch target
            case 'all_conc'
                scores_1 = single_trial_data{comb_bdfs(p,1)}.target{end}.neural_data.dim_red.scores(:,dims);
                fr_2    = single_trial_data{comb_bdfs(p,2)}.target{end}.neural_data.conc_smoothed_fr;
                w       = single_trial_data{comb_bdfs(p,2)}.target{end}.neural_data.dim_red.w(:,dims);
            otherwise
                scores_1 = single_trial_data{comb_bdfs(p,1)}.target{target}.neural_data.dim_red.scores(:,dims);
                fr_2    = single_trial_data{comb_bdfs(p,2)}.target{target}.neural_data.conc_smoothed_fr;          
                w       = single_trial_data{comb_bdfs(p,2)}.target{target}.neural_data.dim_red.w(:,dims);
        end


        % Compute shuffled version of scores for task #2
        % shuffle weights
        w_shuf          = w(randperm(size(w,1)),:);
        scores_2_shuf   = fr_2*w_shuf;
        scores_2_shuf   = scores_2-repmat(mean(scores_2),size(scores_2,1),1);
    
        % do CCA
        [~, ~, t_cc]    = canoncorr( scores_1, scores_2_shuf ); 
        cc(:,i)         = t_cc';
        % and good old corrs
        t_r             = calc_r( scores_1, scores_2_shuf );
        r(:,i)          = sort(abs(t_r),'descend');
  
    end
    
    bootstrap_data(p).cc = cc;
    bootstrap_data(p).r = r;
end




