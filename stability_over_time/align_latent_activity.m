%
% Align latent activity across sessions. Returns a TD object with a new
% field that is the aligned latent activity, and an object with the
% alignment results 
%
%
% ToDo: 
%   - Implement support for multiple areas at the same time
%   - Implement support for Procrustes 
%   - ADDING THE ALIGNED LATENT ACTIVITY IS FOR EACH PAIR OF SESSIONS SO IT
%   CAN'T BE DONE HERE!

function aligned_latent_results = align_latent_activity( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         = [];
xval_yn         = false;
n_folds         = 5;
method          = 'cca'; % 'cca' or 'procrustes'
mani_dims       = 1:10;


if nargin > 1, assignParams(who,params); end % overwrite defaults


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all pairwise comparisons of sessions

sessions        = unique({td.date});
n_sessions      = length(sessions);

% Sometimes the sessions are not sorted by time --fix that here by
% resorting the trials in master_td
[~, i_sort]         = sort(datenum([sessions]));
if sum( i_sort - 1:length(i_sort) ) > 0
   
    sorted_dates = sort( cell2mat( cellfun(@(x) datenum(x), sessions, 'uni', 0) ) );
    for s = 1:n_sessions
        sessions{s} = datestr(sorted_dates(s),'mm-dd-yyyy');
    end
end

comb_sessions   = nchoosek(1:n_sessions,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align using CCA or Procrustes


for c = 1:size(comb_sessions,1)

    
    switch params.xval_yn
    
        % -----------------------------------------------------------------
        % WITHOUT MULTI-FOLD CROSS-VALIDATION
        case false

            % get the data
            [trials1, td1] = getTDidx(td,'date',sessions{comb_sessions(c,1)});
            [trials2, td2] = getTDidx(td,'date',sessions{comb_sessions(c,2)});
            
            
            %trials2 = trials2(randperm(length(trials2)));
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compare dynamics
            aligned_info(c) = compDynamics( td, signals, trials1, trials2, mani_dims, method);
            
            % compare dynamics with good old forrelations
            corr_info(c) = corrDynamics( td, signals, trials1, trials2, mani_dims );
            
            
            % -----------------------------------------------------------------
            % WITH MULTIFOLD CROSS-VALIDATION
        case true
            
            error('Cross-validation not yet implemented when for aligning the dynamics')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Struct to return values 
aligned_latent_results.aligned_info     = aligned_info;
aligned_latent_results.corr_info        = corr_info;
aligned_latent_results.diff_days        = diff_days;
aligned_latent_results.comb_sessions    = comb_sessions;

% add
switch method
    case 'cca'
        aligned_latent_results.cc       = cell2mat(arrayfun(@(x) x.cc, aligned_info, ...
            'uniformoutput', false )');
    case 'procrustes'
        aligned_latent_results.cc       = cell2mat(arrayfun(@(x) x.cc, aligned_info, ...
            'uniformoutput', false )');
end
aligned_latent_results.r                = cell2mat(arrayfun(@(x) x.r, corr_info, ...
                                            'uniformoutput', false )');
aligned_latent_results.signals          = signals;