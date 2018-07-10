%
% Some raw data plots 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot hand kinematics, color-coded per target

% session number
s = 1; %length(sessions);

cols = parula(8);



[~,t_td] = getTDidx(master_td,{'date',sessions{s}});

targets = sort(unique([t_td.target_direction]));
n_targets = numel(targets);

figure
hold on
for t = 1:n_targets
    [~,t_tgt] = getTDidx(t_td,{'target_direction',targets(t)});
    for i = 1:length(t_tgt)
        plot(t_tgt(i).pos(:,1),t_tgt(i).pos(:,2),'color',cols(t,:));
    end
end
set(gca,'TickDir','out','FontSize',14)
xlabel('X position'),ylabel('Y position')
title([master_td(1).monkey ' - ' sessions{s}])

