%
% Removes NaN targets and does other stuff to prepare the data for Han or
% Chips, because of the issues with their TDs
%

if strcmpi(pars.monkey,'han')
    
    % ---------------------------------------------------------------------
    % there are trials with targets = NaN -> exclude them
    master_td       = master_td(~isnan([master_td.target_direction]));
    
    % ---------------------------------------------------------------------
    % in some sessions he had 4 targets and in others. Only keep the common
    % targets
    all_targets     = unique([master_td.target_direction]);
    targets_session = nan(n_sessions,length(all_targets));
    for s = 1:n_sessions
        this_s = getTDidx(master_td,{'date',meta.sessions{s}});
        u_targets = unique([master_td(this_s).target_direction]);
        targets_session(s,1:length(u_targets)) = u_targets;
    end
    
    is_target = zeros(n_sessions,length(all_targets));
    for s = 1:n_sessions
        is_target(s,:) = ismember(all_targets,targets_session(s,:));
    end
    
    % Targets that were present in all sessions
    targets_always_present = all_targets(sum(is_target,1)==n_sessions);
    
    % Only keep trials to targets that were present in all the sessions
    [~,master_td] = getTDidx(master_td,{'target_direction',targets_always_present});

    
    clear unique_targets_session is_target targets_always_present u_targets targets_session this_s all_targets;
% -- And for Chips, to exclude NaN targets
elseif strcmpi(pars.monkey,'chips')
     % there are trials with targets = NaN -> exclude them
    master_td       = master_td(~isnan([master_td.target_direction]));
end