%
% Find optimal trial start for CCA for the reach-to-grasp tasks
%


% PARAMETERS

% trial averaged data?
trial_avg_flg = false;

% bins per trial
n_bins_trial = 50;

% delays we'll try applied on tasks 1 and 2
delays = 1:12;

% manifold dimensionality
mani_dims = 12;

% plot per task
indiv_plot_flg = true;


% ------------------------------------------------------------------------
% load proj params to compute the same CC as in the paper
proj_params = batch_compare_manifold_projs_defaults();

% get bin size
bin_size = datasets{1}.stdata{1}.target{1}.bin_size;

% reach-to-grasp tasks
reach_ds = [4:6 10:11];

% save figs?
save_flg = false;

% where to save figures
dir_figs = '/Users/juangallego/Documents/Publications/2017 - Multi task manifold + some stuff not on DB/additional_analyses_during_review/';

% plot leading mode?
raw_data_plot_flg = true;


%% ------------------------------------------------------------------------
% Do

for d = 1:length(reach_ds)
   
    fname = [datasets{reach_ds(d)}.monkey '_' datasets{reach_ds(d)}.date(1:end-4)];
    disp(fname);
    
    % DO FOR TRIAL-AVERAGED DATA
    if trial_avg_flg 
        
        % The trial averaged data is readily available in the data struct
        sc1 = datasets{reach_ds(d)}.stdata{1}.target{1}.neural_data.dim_red.st_scores_mn(:,1:mani_dims);
        sc2 = datasets{reach_ds(d)}.stdata{2}.target{1}.neural_data.dim_red.st_scores_mn(:,1:mani_dims);
        
        % CC no delay
        win = proj_params.time_win(reach_ds(d),:)/bin_size;
        win_length = win(2)-win(1);
        sc1_nodelay = sc1(win(1):win(2),:);
        sc2_nodelay = sc2(win(1):win(2),:);
        [~,~,cc_no_delay] = canoncorr(sc1_nodelay,sc2_nodelay);
        
        cc_delay1 = zeros(length(delays),mani_dims);
        cc_delay2 = zeros(length(delays),mani_dims);
        
        % do CCA for a range of relative delays
        for dl = 1:length(delays)
            
            % NEED TO PREVENT THIS FROM BREAKING
            
            sc1_delay1 = sc1(delays(dl):delays(dl)+win_length,1:mani_dims);
            sc2_delay1 = sc2(1:win_length+1,1:mani_dims);
            
            if size(sc1_delay1,1)<size(sc2_delay1,1), break; end
            
            [~,~,cc_delay1(dl,:)] = canoncorr(sc1_delay1,sc2_delay1);
        end
        
        for dl = 1:length(delays)
            
            sc1_delay1 = sc1(1:win_length+1,1:mani_dims);
            sc2_delay1 = sc2(delays(dl):delays(dl)+win_length,1:mani_dims);
            
            if size(sc2_delay1,1)<size(sc1_delay1,1), break; end
            
            [~,~,cc_delay2(dl,:)] = canoncorr(sc1_delay1,sc2_delay1);
        end
            
    % DO FOR SINGLE TRIAL DATA
    else

        % Prepare the single trials-- Taken from canon_corr_all_manifolds
        % BUT note that we don't cut them here because we want to play with
        % the delays
        stdata = datasets{reach_ds(d)}.stdata;
        
        % 2) equalize number of trials for all targets of a given task
        for i = 1:length(stdata)
            stdata{i} = equalize_nbr_trials_p_target( stdata{i} );
        end

        % 3) equalize number of trials across tasks
        stdata = equalize_nbr_trials_across_tasks( stdata, 'all_conc' );

        
        % Get the single trial scores in a tensor time x neurons x trials
        ssc1 = stdata{1}.target{end}.neural_data.dim_red.st_scores;
        ssc2 = stdata{2}.target{end}.neural_data.dim_red.st_scores;
        
        
        % -----------------------------------------------------------------
        % CC with no delay in the same window as in the paper
        win = proj_params.time_win(reach_ds(d),:)/bin_size;
        win_length = win(2)-win(1);
        
        
        % get data in this window
        wssc1 = ssc1(win(1):win(2),1:mani_dims,:);
        wssc2 = ssc2(win(1):win(2),1:mani_dims,:);
        
        % turn this into (time x trials) x neurons matrix, taking the time
        % window that we want
        pssc1 = permute(wssc1,[1 3 2]);
        ssc1_nodelay = reshape(pssc1, [], size(pssc1,3));
        
        pssc2 = permute(wssc2,[1 3 2]);
        ssc2_nodelay = reshape(pssc2, [], size(pssc2,3));
        
%         % keep the dynamics within the manifold only
%         ssc1_nodelay = ssc1_nodelay(:,1:mani_dims);
%         ssc2_nodelay = ssc2_nodelay(:,1:mani_dims);
        
        [~,~,cc_no_delay] = canoncorr(ssc1_nodelay,ssc2_nodelay);

        
        % -----------------------------------------------------------------
        % preallocate
        cc_delay1 = nan(length(delays),mani_dims);
        cc_delay2 = nan(length(delays),mani_dims);
        
        
        % do CCA for a range of relative delays
        for dl = 1:length(delays)
            
            if delays(dl)+win_length > min([size(ssc1,1),size(ssc2,1)]), break; end
            
            % get data in this window
            wssc1 = ssc1(delays(dl):delays(dl)+win_length,1:mani_dims,:);
            wssc2 = ssc2(win(1):win(2),1:mani_dims,:);
        
            pssc1 = permute(wssc1,[1 3 2]);
            ssc1_delay1 = reshape(pssc1, [], size(pssc1,3));
            
            pssc2 = permute(wssc2,[1 3 2]);
            ssc2_delay1 = reshape(pssc2, [], size(pssc2,3));
            
            
            if size(ssc1_delay1,1)<size(ssc2_delay1,1), break; end
            
            [~,~,cc_delay1(dl,:)] = canoncorr(ssc1_delay1,ssc2_delay1);
        end
        
        for dl = 1:length(delays)
            
            if delays(dl)+win_length > min([size(ssc1,1),size(ssc2,1)]), break; end
            
            % get data in this window
            wssc1 = ssc1(win(1):win(2),1:mani_dims,:);
            wssc2 = ssc2(delays(dl):delays(dl)+win_length,1:mani_dims,:);
        
            pssc1 = permute(wssc1,[1 3 2]);
            ssc1_delay2 = reshape(pssc1, [], size(pssc1,3));
            
            pssc2 = permute(wssc2,[1 3 2]);
            ssc2_delay2 = reshape(pssc2, [], size(pssc2,3));
            
            
            if size(ssc1_delay2,1)<size(ssc2_delay2,1), break; end
            
            [~,~,cc_delay2(dl,:)] = canoncorr(ssc1_delay2,ssc2_delay2);
        end
    end
    

    % PLOT
    if indiv_plot_flg
        
        cols_plot = parula(numel(delays)+1);
        
        hf = figure('units','normalized','outerposition',[0.3 0.3 0.4 0.4]);
        subplot(121), hold on
        for dl = 1:length(delays)
            plot(cc_delay1(dl,:),'color',cols_plot(dl,:),'linewidth',1)
            lg{dl} = num2str(delays(dl));
        end
        lg{dl+1} = 'no delay';
        plot(cc_no_delay,'k','linewidth',2)
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, 
        ylabel('Canonical correlation'), xlabel('Neural mode')
        title('Delaying modes task 1'), xlim([0 mani_dims]),ylim([0 1])
        if trial_avg_flg
            legend(lg,'Location','SouthWest'),legend boxoff
        else
            legend(lg,'Location','NorthEast'),legend boxoff
        end
        subplot(122), hold on
        for dl = 1:length(delays)
            plot(cc_delay2(dl,:),'color',cols_plot(dl,:),'linewidth',1)
        end
        plot(cc_no_delay,'k','linewidth',2)
        xlabel('Neural mode')
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off, 
        set(hf, 'color', [1 1 1]);
        title('Delaying modes task 2'), xlim([0 mani_dims]),ylim([0 1])
    end
    
    
    % SAVE fig
    if save_flg
        print([dir_figs filesep fname],'-dpng')
    end
    
    
    if raw_data_plot_flg
       
        mn1 = mean(squeeze(ssc1(:,1,:)),2);
        sd1 = std(squeeze(ssc1(:,1,:)),0,2);
        mn2 = mean(squeeze(ssc2(:,1,:)),2);
        sd2 = std(squeeze(ssc2(:,1,:)),0,2);
        
        hrf = figure('units','normalized','outerposition',[0.3 0.3 0.4 0.4]);
        subplot(121), hold on
        for i = 1:size(ssc1,3), plot(squeeze(ssc1(:,1,i)),'c'); end
        plot(mn1,'b','linewidth',2)
        plot(mn1+sd1,'-.b','linewidth',2)
        plot(mn1-sd1,'-.b','linewidth',2)
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
        xlabel('Bin nbr'),ylabel('Neural mode 1')
        title(datasets{reach_ds(d)}.labels{1})
        subplot(122),hold on
        for i = 1:size(ssc2,3), plot(squeeze(ssc2(:,1,i)),'color',[.65 .65 .65]); end
        plot(mn2,'k','linewidth',2)
        plot(mn2+sd2,'-.k','linewidth',2)
        plot(mn2-sd2,'-.k','linewidth',2)
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
        xlabel('Bin nbr'),ylabel('Neural mode 2')
        title(datasets{reach_ds(d)}.labels{2})
        set(hrf, 'color', [1 1 1]);
    end
end