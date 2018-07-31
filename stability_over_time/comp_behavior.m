%
% Compare kinematics
%

function results = comp_behavior( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signal                  = [];
n_folds                 = 6; % 'cca' or 'procrustes'
trial_avg               = false;


if nargin > 1, assignParams(who,params); end % overwrite defaults


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
n_comb_sessions         = size(comb_sessions,1);

n_signals               = size(td(1).(signal),2);

if trial_avg
    
    td                  = trialAverage( td, {'target_direction','date'} );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

r                       = zeros(n_comb_sessions,n_signals);

disp('Comparing similarity kinematics');

for c = 1:n_comb_sessions

    % get two TD structs, one per session to compare
    [trials1, ~]        = getTDidx(td,'date',sessions{comb_sessions(c,1)});
    [trials2, ~]        = getTDidx(td,'date',sessions{comb_sessions(c,2)});

    
    tr                  = corrDynamics( td, signal, trials1, trials2, 1:n_signals );
    r(c,:)              = tr.r;
end


% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end
results.diff_days           = diff_days;
results.comb_sessions       = comb_sessions;
results.r                   = r;
results.var                 = signal;