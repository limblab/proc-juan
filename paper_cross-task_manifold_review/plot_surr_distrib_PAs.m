%
% Plot surrogate PA distributions -> Need to generate the data first with
% signif_PAs_using_TME.m
%


% Dataset number
ds_example = 1;

% What do you want to do?
plot_example = false;
plot_all_dims = true;


% PLOT ONE EXAMPLE
if plot_example
    
    t_surrs = all_TME_PAs{ds_example};

    n_PAs = size(t_surrs,2);
    n_surrs = size(t_surrs,1);

    % compute distributions
    x_hist = 0:pi/2/900:pi/2;
    y_hist = zeros(length(x_hist)-1,n_PAs );

    for i = 1:n_PAs 
        y_hist(:,i) = histcounts(t_surrs(:,i),x_hist) / n_surrs * 100;
    end



    % Plot
    cols = hot(n_PAs+5);

    figure, hold on
    for i = 1:n_PAs   
        plot(rad2deg(x_hist(1:end-1)),y_hist(:,i),'color',cols(i,:),'linewidth',1)
    end
    set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(gcf, 'color', [1 1 1]);
    xlabel('Principal angle (deg)'), ylabel('Manifold comparisons (%)');
end



% PLOT THRESHOLD FOR ALL DATASETS
if plot_all_dims
   
    % number of sessions (=nbr different space dimensionalities)
    n_sessions = length(datasets);
    % manifold dimensionality
    n_PAs = size(TME_th,2);
    % dimensionality of the neural space for each session:
    dimens_p_session = cellfun(@(x) length(x.neural_chs), datasets);
    
    % to save the mean TME threshold per space dimensionality
    mean_TME_th_p_dimens = zeros(n_sessions,n_PAs);
%     sd_TME_th_p_dimens = zeros(n_sessions,n_PAs);
    
    % fill mean_TME_th_p_dimens
    for d = 1:length(datasets)
        
        idx_t_session = find(session_nbr == d)
        mean_TME_th_p_dimens(d,:) = mean(TME_th(idx_t_session,:),1);
%         sd_TME_th_p_dimens(d,:) = std(TME_th(idx_t_session,:),0,1);
    end
    
    % Find and sort unique space dimensionalities
    [sorted_dimens, idx_sorted_dimens] = unique(dimens_p_session);
    
    % Resort threshold based on the dimensionalities of the neural space
    sorted_mean_TME_th_p_dimens = mean_TME_th_p_dimens(idx_sorted_dimens,:);
    
    unique_dimens = length(sorted_dimens);
    
    % Plot
    cols = parula(unique_dimens +2);
    for i = 1:unique_dimens
       lgnd_dim{i} = ['N=' num2str(sorted_dimens(i))];
    end
    
    figure,hold on
    for i = 1:unique_dimens 
        plot(rad2deg(sorted_mean_TME_th_p_dimens(i,:)),'color',cols(i,:),'linewidth',1)
    end
    set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, set(gcf, 'color', [1 1 1]);
    xlabel('Neural mode'), ylabel(['Princ. angle for P<' num2str(P_orth)]);
    ylim([0 90]),legend(lgnd_dim,'location','southeast'),legend boxoff
    
    
end