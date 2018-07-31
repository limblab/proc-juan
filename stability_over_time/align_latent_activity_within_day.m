%
%
%

function within_day_results = align_latent_activity_within_day( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         = [];
xval_yn         = false;
n_shuff         = 100;
method          = 'cca'; % 'cca' or 'procrustes'
mani_dims       = 1:10;


if nargin > 1, assignParams(who,params); end % overwrite defaults


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get sessions

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align


for s = 1:n_sessions
    
    % get all trials in this session
    [~, tdd]    = getTDidx(td,'date',sessions{s});

    % get targets in this session, and how many times they appear
    tgt         = [tdd.target_direction];
    tgts        = unique(tgt);
    n_tgts      = length(unique(tgt));
    
    reps_p_tgt  = zeros(1,n_tgts);
    for t = 1:n_tgts
       reps_p_tgt(t) = sum(tgt == tgts(t));
    end
    
    % we're going to divide the trials in two blocks, see how many to keep
    % (min of trials for all the targets)
    trials_p_block = floor( min(reps_p_tgt) / 2 );
    
    
    % Compare randomly selected subsets of trials
    for r = 1:n_shuff
        
        % randomize trial order and define two non-overlapping subsets of
        % trials
        idx_rand    = randperm(length(tgt));
        
        rnd_tgt     = tgt(idx_rand);
        
        idx_ss1     = zeros(1, trials_p_block * n_tgts );
        idx_ss2     = zeros(1, trials_p_block * n_tgts );
        
        for t = 1:n_tgts
           
            idx_t_tgt = find( rnd_tgt == tgts(t) );
            
            unshuf_idx = idx_rand( idx_t_tgt ); % this puts the indices back in the real trial "coordinate system"
            
            idx_ss1( (t-1)*trials_p_block+1 : t*trials_p_block ) = ...
                unshuf_idx( 1 : trials_p_block );
            idx_ss2( (t-1)*trials_p_block+1 : t*trials_p_block ) = ...
                unshuf_idx( end-trials_p_block+1 : end );
        end
        
        aligned_info(r) = compDynamics( tdd, signals, idx_ss1, idx_ss2, mani_dims );
    end
    
    
    % Summarize results
    pooled_ccs      = cell2mat( arrayfun(@(x) x.cc', aligned_info,'uni',false) )';
    
    
    % Save
    within_day_results(s).aligned_info = aligned_info;
    within_day_results(s).session_idx = s;
    within_day_results(s).session = sessions{s};
    within_day_results(s).cc_m = mean(pooled_ccs,1);
    within_day_results(s).cc_sd = std(pooled_ccs,0,1);
    within_day_results(s).cc_prc99 = prctile(pooled_ccs,99);
end


% %% -----------------------------------------------------------------------------------------------------
% % PLOT
% 
% top_lv_plot = 4;
% 
% ccs_plot = cell2mat( arrayfun( @(x) x.cc_prc99, within_day_align_results, 'uni', false )' );
% 
% m_ccs_plot = mean(ccs_plot(:,1:top_lv_plot),2);
% sem_ccs_plot = mean(ccs_plot(:,1:top_lv_plot),2)/sqrt(size(ccs_plot,1));
% 
% 
% figure,errorbar(m_ccs_plot',sem_ccs_plot','.k','markersize',32,'linestyle','none')
% ylim([0 1]),xlim([0 n_sessions+1])
% set(gca,'TickDir','out','FontSize',14), box off