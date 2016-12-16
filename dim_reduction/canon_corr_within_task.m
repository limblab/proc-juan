%
% Canonical correlation among trials of the same session, to assess
% inter-trial similarity. Also computes bootstrapped confidence bounds
% (99th percentile of 1,000 shuffles)
%
%   cc_within = canon_corr_within_task( single_trial_data, dims, varargin )
%
%
% Inputs(opt)           : [default]
%   single_trial_data   : single_trial_data struct for one session
%   dims                : manifold dimensions --a vector, e.g. 1:10
%   (shuff)             : [false] shuffle trial order 
%   (target)            : ['all_conc'] target or targets that will be
%                           analysed
%   (time_window)       : [0 0.5] portion of the trial that will be used fo
%                           the analysis in s
%   (conf_thr)          : [true] get bootstrapped confidence threshold
%
% Outputs:
%   cc_within           : struct with CC, linear transformations, projected
%                           latent variables and bootstrapped conf. bounds
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
if nargin >= 5
    time_window     = varargin{3};
    if size(time_window) ~= [1, 2], error('time_window has wrong dimensions'); end
else
    time_window     = [0 0.5];
end
if nargin == 6
    conf_thr        = varargin{4};
else
    conf_thr        = true;
end

                        
% % -------------------------------------------------------------------------
% % get some meta info
% nbr_bdfs            = length(single_trial_data);
% 
% 
% % -------------------------------------------------------------------------
% % prepare data for doing canonical correlation
% 
% % 1) equalize trial duration across all tasks
% single_trial_data   = equalize_single_trial_dur( single_trial_data, ...
%                         'time_win', time_window );
                        

% -------------------------------------------------------------------------
% do

   
split_data          = split_single_trial_data( single_trial_data, 2, shuff );

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


% -------------------------------------------------------------------------
% if want to obtain the significance threshold with bootstrapping
if conf_thr
    % do bootstrapping with default params -1,000 repetitions and 99th
    % percentile
    signif_boots                = bootstrap_canon_corr( scores_1, scores_2 );
end


% -------------------------------------------------------------------------
% store, to return
cc_within.lin_transform.A       = A;
cc_within.lin_transform.B       = B;
cc_within.cc                    = r;
cc_within.lin_transform.U       = U;
cc_within.lin_transform.V       = V;
cc_within.stats                 = stats;
if exist('signif_boots','var')
    cc_within.signif_boots      = signif_boots;
end