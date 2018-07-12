%
% Compare kinematics
%

function results = comp_behavior( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signal                  = [];
n_folds                 = 6; % 'cca' or 'procrustes'


if nargin > 1, assignParams(who,params); end % overwrite defaults


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STUFF

% get all pairs of sessions
sessions                = unique({td.date});
n_sessions              = length(sessions);
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);

n_signals               = size(td(1).(signal),2);

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