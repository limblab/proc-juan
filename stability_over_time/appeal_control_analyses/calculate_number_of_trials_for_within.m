


% get the number of trials that will be used within-day
dates = unique({master_td_all_trials.date});

n_trials = zeros(1,length(dates));
for iDate = 1:length(dates)
    [~,td] = getTDidx(master_td_all_trials,'date',dates{iDate});
    
        % get targets in this session, and how many times they appear
    tgt         = [td.target_direction];
    tgts        = unique(tgt);
    n_tgts      = length(unique(tgt));
    
    reps_p_tgt  = zeros(1,n_tgts);
    for t = 1:n_tgts
       reps_p_tgt(t) = sum(tgt == tgts(t));
    end
    
    n_trials(iDate) = floor( min(reps_p_tgt) / 2 );
    
end



