%
% Plot neural variance per marginalization, pooled over all sessions for
% each monkey separately
%
% Updated during the review, to show individual datapoints
%


% datasets per monkey
s_p_monkey      = {1:3,4:6}; % C and J

% paper colors
marg_cols       = [.4 .6 .5; .7 .1 .1; .5 .65 .9; 1 .6 .3];

% marg labels
marg_labels     = {'task','target','time','task-target'};

% color individual datapoints
pt_cols         = [.75 .75 .75];


% create matrices with data
for m = 1:length(s_p_monkey)
    
    % cell array expl_var per monkey (session x marg )
    expl_var{m} = zeros(length(s_p_monkey{m}),length(marg_labels));
    
    for s = 1:length(s_p_monkey{m})
        
        expl_var{m}(s,:) = dPCA_results{s_p_monkey{m}(s)}.expl_var.totalMarginalizedVar / ...
                        dPCA_results{s_p_monkey{m}(s)}.expl_var.totalVar * 100;
    end
end


% do stats
m_var_marg  = cell2mat(cellfun(@(x) mean(x), expl_var,'UniformOutput',false)');
sd_var_marg = cell2mat(cellfun(@(x) std(x), expl_var,'UniformOutput',false)');


% change order of the data and the lab so marginalizations are shown as
% task-independet and task dependent
marg_swap   = [3 2 1 4];
m_var_marg  = m_var_marg(:,marg_swap);
sd_var_marg = sd_var_marg(:,marg_swap);
marg_labels = marg_labels(marg_swap);
marg_cols   = marg_cols(marg_swap,:);
% rearrange order individual datapoints according to marg_swap
for k = 1:length(s_p_monkey)
    expl_var{k} = expl_var{k}(:,marg_swap);
end


% plot
figure, hold on
for m = 1:length(marg_labels)
    for k = length(s_p_monkey):-1:1
        if k == 1 % for monkey 1 empty bars
            errorbar( 3*m-1, m_var_marg(k,m), 0, sd_var_marg(k,m),'color',marg_cols(m,:),'marker','none','Linewidth',1.5);
            tb(m) = bar( 3*m-1, m_var_marg(k,m) );
            set(tb(m),'FaceColor',[1 1 1],'EdgeColor',marg_cols(m,:),'LineWidth',3)
            % plot individual datapoints
%            plot(3*m-1+.25, expl_var{k}(:,m), 'o', 'color', marg_cols(m,:), 'linestyle', 'none')
            plot(3*m-1+.25, expl_var{k}(:,m), 'o', 'color', pt_cols, 'linestyle', 'none')
        else % for monkey 2 filled bars
            errorbar( 3*m-2, m_var_marg(k,m), 0, sd_var_marg(k,m),'color',marg_cols(m,:),'marker','none','Linewidth',1.5);
            tb2(m) = bar( 3*m-2, m_var_marg(k,m) );
            set(tb2(m),'FaceColor',marg_cols(m,:),'EdgeColor',marg_cols(m,:),'LineWidth',3)
            % plot individual datapoints
%             plot(3*m-2+.25, expl_var{k}(:,m), 'o', 'color', marg_cols(m,:), 'linestyle', 'none')
            plot(3*m-2+.25, expl_var{k}(:,m), 'o', 'color', pt_cols, 'linestyle', 'none')
        end
    end
end
set(gca,'TickDir','out','FontSize',14), box off
ylabel('Neural var. expl. (%)'), ylim([0 45])
set(gca,'XTick',1.5:3:10.5,'XTickLabel',marg_labels,'XTickLabelRotation',45)
legend([tb2(1) tb(1)],'J','C'),legend boxoff


% Clear vars
clearvars -except *_results *_params datasets manifold_dim