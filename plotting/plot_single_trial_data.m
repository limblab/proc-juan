%
% Plot single trial data struct
%

function plot_single_trial_data(single_trial_data )




% make these parameters
neural_PCs              = 1:6;
neural_PCs_vs_muscles   = 1:3;
muscle_PCs              = 1:min(3,size(dim_red_emg.scores,2));

plot_all                = false; % plot all or just mean +/- SD


colors                  = parula(nbr_targets);

% neural PCs only
nbr_rows                = floor(sqrt(length(neural_PCs)));
nbr_cols                = ceil(length(neural_PCs)/nbr_rows);
max_score               = max(cell2mat(cellfun(@(x) max(max(max(abs(x.neural_scores.data(:,neural_PCs,:))))),...
                            single_trial_data,'UniformOutput',false))); % improve coding !!

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(neural_PCs)
   subplot(nbr_rows,nbr_cols,i), hold on
   for t = 1:nbr_targets
        t_axis           = 0:bin_size:(size(single_trial_data{t}.neural_scores.data,1)...
                            -1)*bin_size;
        if plot_all
            plot(t_axis,squeeze(single_trial_data{t}.neural_scores.data(:,neural_PCs(i),:)),...
                'color',colors(t,:));
        else
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs(i)),...
                            'color',colors(t,:),'linewidth',6);
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs(i)) + ...
                single_trial_data{t}.neural_scores.sd(:,neural_PCs(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs(i)) - ...
                single_trial_data{t}.neural_scores.sd(:,neural_PCs(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
        end
        clear t_axis;
        ylim([-max_score, max_score])
   end
   set(gca,'TickDir','out'),set(gca,'FontSize',14); 
   title(['neural comp ' num2str(neural_PCs(i))]);
   if i == 1 || rem(i-1,nbr_cols) == 0
       ylabel('a.u.','FontSize',14)
   end
   if i > nbr_cols*(nbr_rows-1)
       xlabel('time (s)','FontSize',14)
   end
end



% neural PCs vs muscles
nbr_rows                = floor(sqrt(length([neural_PCs_vs_muscles, muscle_PCs])));
nbr_cols                = ceil(length([neural_PCs_vs_muscles, muscle_PCs])/nbr_rows);

max_score               = max(cell2mat(cellfun(@(x) max(max(max(abs(x.neural_scores.data(:,neural_PCs,:))))),...
                            single_trial_data,'UniformOutput',false))); % improve coding !!

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(neural_PCs_vs_muscles)
   subplot(nbr_rows,nbr_cols,i), hold on
   for t = 1:nbr_targets
        t_axis           = 0:bin_size:(size(single_trial_data{t}.neural_scores.data,1)...
                            -1)*bin_size;
        if plot_all
            plot(t_axis,squeeze(single_trial_data{t}.neural_scores.data(:,neural_PCs_vs_muscles(i),:)),...
                'color',colors(t,:));
        else
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)),...
                            'color',colors(t,:),'linewidth',6);
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) + ...
                single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
            plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) - ...
                single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
        end
        clear t_axis;
        ylim([-max_score, max_score])
   end
   set(gca,'TickDir','out'),set(gca,'FontSize',14); 
   title(['neural comp ' num2str(neural_PCs_vs_muscles(i))]);
   if i == 1 || rem(i-1,nbr_cols) == 0
       ylabel('a.u.','FontSize',14)
   end
end

for i = 1:length(muscle_PCs)
   subplot(nbr_rows,nbr_cols,i+length(neural_PCs_vs_muscles)), hold on
   for t = 1:nbr_targets
        t_axis           = 0:bin_size:(size(single_trial_data{t}.emg_scores.data,1)...
                            -1)*bin_size;
        if plot_all
            plot(t_axis,squeeze(single_trial_data{t}.emg_scores.data(:,muscle_PCs(i),:)),...
                'color',colors(t,:));
        else
            plot(t_axis,single_trial_data{t}.emg_scores.mn(:,muscle_PCs(i)),...
                            'color',colors(t,:),'linewidth',6);
            plot(t_axis,single_trial_data{t}.emg_scores.mn(:,muscle_PCs(i)) + ...
                single_trial_data{t}.emg_scores.sd(:,muscle_PCs(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
            plot(t_axis,single_trial_data{t}.emg_scores.mn(:,muscle_PCs(i)) - ...
                single_trial_data{t}.emg_scores.sd(:,muscle_PCs(i)), ...
                'color',colors(t,:),'linewidth',2,'linestyle','-.');
        end
        clear t_axis;
   end
   set(gca,'TickDir','out'),set(gca,'FontSize',14); 
   title(['muscle comp ' num2str(muscle_PCs(i))]);
   if i == 1 || rem(i-1,nbr_cols) == 0
       ylabel('a.u.','FontSize',14)
   end
   xlabel('time (s)','FontSize',14)
end




% neural PCs vs position
if ~isempty(cropped_binned_data.cursorposbin)
                    
    nbr_rows                = floor(sqrt(length([neural_PCs_vs_muscles, muscle_PCs])));
    nbr_cols                = ceil(length([neural_PCs_vs_muscles, muscle_PCs])/nbr_rows);

    max_score               = max(cell2mat(cellfun(@(x) max(max(max(abs(x.neural_scores.data(:,neural_PCs,:))))),...
                                single_trial_data,'UniformOutput',false))); % improve coding !!

    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:length(neural_PCs_vs_muscles)
       subplot(nbr_rows,nbr_cols,i), hold on
       for t = 1:nbr_targets
            t_axis           = 0:bin_size:(size(single_trial_data{t}.neural_scores.data,1)...
                                -1)*bin_size;
            if plot_all
                plot(t_axis,squeeze(single_trial_data{t}.neural_scores.data(:,neural_PCs_vs_muscles(i),:)),...
                    'color',colors(t,:));
            else
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)),...
                                'color',colors(t,:),'linewidth',6);
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) + ...
                    single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) - ...
                    single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
            end
            clear t_axis;
            ylim([-max_score, max_score])
       end
       set(gca,'TickDir','out'),set(gca,'FontSize',14); 
       title(['neural comp ' num2str(neural_PCs_vs_muscles(i))]);
       if i == 1 || rem(i-1,nbr_cols) == 0
           ylabel('a.u.','FontSize',14)
       end
    end

    for i = 1:size(single_trial_data{1}.pos.data,2)
       subplot(nbr_rows,nbr_cols,i+length(neural_PCs_vs_muscles)), hold on
       for t = 1:nbr_targets
            t_axis           = 0:bin_size:(size(single_trial_data{t}.pos.data,1)...
                                -1)*bin_size;
            if plot_all
                plot(t_axis,squeeze(single_trial_data{t}.pos.data(:,i,:)),...
                    'color',colors(t,:));
            else
                plot(t_axis,single_trial_data{t}.pos.mn(:,i),...
                                'color',colors(t,:),'linewidth',6);
                plot(t_axis,single_trial_data{t}.pos.mn(:,i) + ...
                    single_trial_data{t}.pos.sd(:,muscle_PCs(i)), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
                plot(t_axis,single_trial_data{t}.pos.mn(:,i) - ...
                    single_trial_data{t}.pos.sd(:,i), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
            end
            clear t_axis;
       end
       set(gca,'TickDir','out'),set(gca,'FontSize',14); 
       title(['pos/force ' num2str(i)]);
       if i == 1 || rem(i-1,nbr_cols) == 0
           ylabel('a.u.','FontSize',14)
       end
       xlabel('time (s)','FontSize',14)
    end


    % neural PCs vs velocity
    nbr_rows                = floor(sqrt(length([neural_PCs_vs_muscles, muscle_PCs])));
    nbr_cols                = ceil(length([neural_PCs_vs_muscles, muscle_PCs])/nbr_rows);

    max_score               = max(cell2mat(cellfun(@(x) max(max(max(abs(x.neural_scores.data(:,neural_PCs,:))))),...
                                single_trial_data,'UniformOutput',false))); % improve coding !!

    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:length(neural_PCs_vs_muscles)
       subplot(nbr_rows,nbr_cols,i), hold on
       for t = 1:nbr_targets
            t_axis           = 0:bin_size:(size(single_trial_data{t}.neural_scores.data,1)...
                                -1)*bin_size;
            if plot_all
                plot(t_axis,squeeze(single_trial_data{t}.neural_scores.data(:,neural_PCs_vs_muscles(i),:)),...
                    'color',colors(t,:));
            else
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)),...
                                'color',colors(t,:),'linewidth',6);
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) + ...
                    single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
                plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) - ...
                    single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
            end
            clear t_axis;
            ylim([-max_score, max_score])
       end
       set(gca,'TickDir','out'),set(gca,'FontSize',14); 
       title(['neural comp ' num2str(neural_PCs_vs_muscles(i))]);
       if i == 1 || rem(i-1,nbr_cols) == 0
           ylabel('a.u.','FontSize',14)
       end
    end

    for i = 1:size(single_trial_data{1}.vel.data,2)
       subplot(nbr_rows,nbr_cols,i+length(neural_PCs_vs_muscles)), hold on
       for t = 1:nbr_targets
            t_axis           = 0:bin_size:(size(single_trial_data{t}.vel.data,1)...
                                -1)*bin_size;
            if plot_all
                plot(t_axis,squeeze(single_trial_data{t}.vel.data(:,i,:)),...
                    'color',colors(t,:));
            else
                plot(t_axis,single_trial_data{t}.vel.mn(:,i),...
                                'color',colors(t,:),'linewidth',6);
                plot(t_axis,single_trial_data{t}.vel.mn(:,i) + ...
                    single_trial_data{t}.vel.sd(:,i), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
                plot(t_axis,single_trial_data{t}.vel.mn(:,i) - ...
                    single_trial_data{t}.vel.sd(:,i), ...
                    'color',colors(t,:),'linewidth',2,'linestyle','-.');
            end
            clear t_axis;
       end
       set(gca,'TickDir','out'),set(gca,'FontSize',14); 
       title(['velocity ' num2str(i)]);
       if i == 1 || rem(i-1,nbr_cols) == 0
           ylabel('a.u.','FontSize',14)
       end
       xlabel('time (s)','FontSize',14)
    end

end

% % neural PCs vs some EMGs
% emgs_to_plot            = 10:12;
% nbr_rows                = floor(sqrt(length([neural_PCs_vs_muscles, emgs_to_plot])));
% nbr_cols                = ceil(length([neural_PCs_vs_muscles, emgs_to_plot])/nbr_rows);
% 
% max_score               = max(cell2mat(cellfun(@(x) max(max(max(abs(x.neural_scores.data(:,neural_PCs,:))))),...
%                             single_trial_data,'UniformOutput',false))); % improve coding !!
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for i = 1:length(neural_PCs_vs_muscles)
%    subplot(nbr_rows,nbr_cols,i), hold on
%    for t = 1:nbr_targets
%         t_axis           = 0:bin_size:(size(single_trial_data{t}.neural_scores.data,1)...
%                             -1)*bin_size;
%         if plot_all
%             plot(t_axis,squeeze(single_trial_data{t}.neural_scores.data(:,neural_PCs_vs_muscles(i),:)),...
%                 'color',colors(t,:));
%         else
%             plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)),...
%                             'color',colors(t,:),'linewidth',6);
%             plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) + ...
%                 single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
%                 'color',colors(t,:),'linewidth',2,'linestyle','-.');
%             plot(t_axis,single_trial_data{t}.neural_scores.mn(:,neural_PCs_vs_muscles(i)) - ...
%                 single_trial_data{t}.neural_scores.sd(:,neural_PCs_vs_muscles(i)), ...
%                 'color',colors(t,:),'linewidth',2,'linestyle','-.');
%         end
%         clear t_axis;
%         ylim([-max_score, max_score])
%    end
%    set(gca,'TickDir','out'),set(gca,'FontSize',14); 
%    title(['neural comp ' num2str(neural_PCs_vs_muscles(i))]);
%    if i == 1 || rem(i-1,nbr_cols) == 0
%        ylabel('a.u.','FontSize',14)
%    end
% end
% 
% for i = 1:numel(emgs_to_plot)
%    subplot(nbr_rows,nbr_cols,i+length(neural_PCs_vs_muscles)), hold on
%    for t = 1:nbr_targets
%         t_axis           = 0:bin_size:(size(single_trial_data{t}.vel.data,1)...
%                             -1)*bin_size;
%         if plot_all
%             plot(t_axis,squeeze(single_trial_data{t}.vel.data(:,i,:)),...
%                 'color',colors(t,:));
%         else
%             plot(t_axis,single_trial_data{t}.vel.mn(:,i),...
%                             'color',colors(t,:),'linewidth',6);
%             plot(t_axis,single_trial_data{t}.vel.mn(:,i) + ...
%                 single_trial_data{t}.vel.sd(:,muscle_PCs(i)), ...
%                 'color',colors(t,:),'linewidth',2,'linestyle','-.');
%             plot(t_axis,single_trial_data{t}.vel.mn(:,i) - ...
%                 single_trial_data{t}.vel.sd(:,i), ...
%                 'color',colors(t,:),'linewidth',2,'linestyle','-.');
%         end
%         clear t_axis;
%    end
%    set(gca,'TickDir','out'),set(gca,'FontSize',14); 
%    title(['velocity ' num2str(i)]);
%    if i == 1 || rem(i-1,nbr_cols) == 0
%        ylabel('a.u.','FontSize',14)
%    end
%    xlabel('time (s)','FontSize',14)
% end