%
% Do output potent/null analysis of all tasks in a session. 
%
%
%

function opn_spaces = output_potent_null_spaces_all_manifolds( single_trial_data, varargin )


% read input params, if passed
if nargin == 2
    params              = output_potent_null_spaces_defaults( varargin{1} );
else
    params              = output_potent_null_spaces_defaults;
end


nbr_bdfs                = length(single_trial_data);

% -------------------------------------------------------------------------
% prepare data for the analysis

% % 1) equalize trial duration across all tasks
% single_trial_data       = equalize_single_trial_dur( single_trial_data, ...
%                             'time_win', params.time_win );
% 
% % 2) equalize number of trials for all targets of a given task
% for i = 1:nbr_bdfs
%     single_trial_data{i} = equalize_nbr_trials_p_target( single_trial_data{i} );
% end
% 
% % 3) equalize number of trials across tasks
% single_trial_data       = equalize_nbr_trials_across_tasks( single_trial_data, params.target );

 
% -------------------------------------------------------------------------
% do !

for i = 1:nbr_bdfs
    opn_spaces(i)       = output_potent_null_spaces( single_trial_data{i}, params );
end


% -------------------------------------------------------------------------
% plots

if params.detailed_plots_yn
    
    nbr_rows            = floor(sqrt(nbr_bdfs));
    nbr_cols            = ceil(nbr_bdfs/nbr_rows);
    nbr_emgs            = length(single_trial_data{1}.target{1}.emg_data.emg_names);
    
    colors              = distinguishable_colors(nbr_emgs);
    figure
    for i = 1:nbr_bdfs
        
        mean_R2         = mean(opn_spaces(i).stats_W(:,2:end),2);
        std_R2          = std(opn_spaces(i).stats_W(:,2:end),0,2);
        
        subplot(nbr_rows,nbr_cols,i), hold on
        for e = 1:nbr_emgs
            % plot not-crossvalidated R^2
            plot(0,opn_spaces(i).stats_W(e,1),'v','markersize',14,'color',colors(e,:),'linewidth',3)
            % and mean +/-SD of the cross-validated EMGs
            errorbar(e,mean_R2(e),std_R2(e),'marker','o','markersize',14,'color',colors(e,:),'linewidth',3)
        end
        set(gca,'TickDir','out','FontSize',14)
        xlim([-1 nbr_emgs + 1]),ylim([-0.1 1.1])
        set(gca,'XTick',1:nbr_emgs), set(gca,'XTickLabel',single_trial_data{1}.target{1}.emg_data.emg_names)
        set(gca,'XTickLabelRotation',45)
        ylabel('R^2 EMG predictions')
    end
end