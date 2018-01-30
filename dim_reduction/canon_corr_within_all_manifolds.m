%
%
%


function cc_within = canon_corr_within_all_manifolds( single_trial_data, ...
                        dims, varargin )

                        
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
    nbr_shuff       = varargin{4};
else
    % the default nbr of shuffles depends on whether we want to compare two
    % halves of a task or two sets of random trials
    if shuff
        nbr_shuff   = 100;
    else 
        nbr_shuff   = 1;
    end
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

for i = 1:nbr_bdfs
   
    % if comparing random subsets, do n times, otherwise do 1 time
    for s = 1:nbr_shuff       
        
        cc_data{i}.data{s} = canon_corr_within_task( single_trial_data{i}, ...
                        dims, shuff, target, time_window );
    end
end


% -------------------------------------------------------------------------
% add some stats (mean, SD, 99th percentile)

if shuff
    for i = 1:nbr_bdfs
        if nbr_shuff > 1
            cc_within.task{i}.data   = cell2mat( cellfun( @(x) x.cc, cc_data{i}.data, ...
                                    'UniformOutput', false )' );
            cc_within.task{i}.mn     = mean(cc_within.task{i}.data);
            cc_within.task{i}.sd     = std(cc_within.task{i}.data);
            cc_within.task{i}.pctile99 = prctile(cc_within.task{i}.data,99);
        else
            cc_within.task{i}.data   = cc_data{i}.data{1};
            cc_within.task{i}.mn     = mean(cc_data{i}.data{1});
            cc_within.task{i}.sd     = std(cc_data{i}.data{1});
            cc_within.task{i}.pctile99 = prctile(cc_data{i}.data{1},99);
        end
    end
end

% -------------------------------------------------------------------------
% summary stats per pair

bdf_pairs           = nchoosek(1:nbr_bdfs,2);

for p = 1:size(bdf_pairs,1)
    cc_within.pair{p}.mn        = mean([cc_within.task{bdf_pairs(p,1)}.mn; ...
                                    cc_within.task{bdf_pairs(p,2)}.mn]);
    cc_within.pair{p}.pctile99  = mean([cc_within.task{bdf_pairs(p,1)}.pctile99; ...
                                    cc_within.task{bdf_pairs(p,2)}.pctile99]);
end
