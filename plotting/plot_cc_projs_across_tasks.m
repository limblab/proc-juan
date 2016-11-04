%
% Plot CCs between two tasks
%

function plot_cc_projs_across_tasks( proj_results, params, session, pair, varargin )


if nargin >= 5
    nbr_projs       = varargin{1};
else
    nbr_projs       = 3;
end

if nargin >= 6
    bin_size        = varargin{2};
else
    bin_size        = 0.02;
end


% create time vector
t                   = 0:bin_size:bin_size*(length(proj_results.data{session}.can_corrs.lin_transform(1).U)-1);

% define time for vertical lines to split trials
trial_length        = (params.time_win(session,2)-params.time_win(session,1))/bin_size+1;
nbr_trials          = length(t)/trial_length;
vline_samples       = (0:nbr_trials)*trial_length;
vline_samples(1)    = [];

% get axis scale
max_scale           = max(max([abs(proj_results.data{session}.can_corrs.lin_transform(pair).U(:,1:nbr_projs)), ...
                        abs(proj_results.data{session}.can_corrs.lin_transform(pair).V(:,1:nbr_projs))]));
max_scale           = ceil(max_scale);

figure
for i = 1:nbr_projs
    subplot(nbr_projs,1,i), hold on
    plot( t, proj_results.data{session}.can_corrs.lin_transform(pair).U(:,i),'b','linewidth',2)
    plot( t, proj_results.data{session}.can_corrs.lin_transform(pair).V(:,i),'c','linewidth',2)
    for l = 1:length(vline_samples)
        plot([t(round(vline_samples(l))) t(round(vline_samples(l)))],[-max_scale max_scale],'w','linewidth',4)
    end
    
    legend(['r=' num2str(proj_results.data{session}.can_corrs.cc(pair,i))]), legend boxoff
    set(gca,'TickDir','out','FontSize',14)
    xlim([15 30]), ylim([-max_scale, max_scale])
    if i == nbr_projs
        xlabel('time (s)')
    end
    ylabel(['proj. ' num2str(i)])
end
