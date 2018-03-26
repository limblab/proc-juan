%
% Do CCA after aligning the single trials to movement onset, identified
% from the EMGs
%


clearvars -except datasets; clc; close all;


% ------------------------------------------------------------------------
% Main Params

% bins per trial
n_bins_trial = 40; % 50 works well but excludes too many trials

% trim original data window so it has the same length as the new one?
trim_original_win_flg = true;

% manifold dimensionality
mani_dims = 12;

% do control?
do_control_flg = false;
nbr_shuffles = 1000; % nbr shuffles for control
P_th = 0.01;


% ------------------------------------------------------------------------
% Extra plotting

% plot single trials?
plot_single_trial_flg = false;
% pause for each single trial plot (in s)
sec_pause = 0.1;

% plot per session
plot_each_session_flg = true;


% ------------------------------------------------------------------------
% Params for the analysis

% To exclude trials based on weird movement onsets
min_nbr_bins_beginning_mov_onset = 15; % 15 works well
min_nbr_bins_mov_onset_end = 30;
max_diff_mov_onset = 30;

% where to start a trial
trial_begins = 10; % min_nbr_bins_beginning_mov_onset - 1;

% Params for EMG based detection of movement onset
max_emg_beginning_trial = 0.4;
snr_th_emg = 0.45;
fc_lp_emg = 3;
max_emg_mov_onset = 0.10;


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
dir_figs = '/Users/juangallego/Documents/Publications/2017 - Multi task manifold + some stuff not on DB/additional_analyses_during_review/CCA_EMG_dynamics_onset/';




%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% FIND MOVEMENT ONSETS


% preallocate variables
all_r = [];
all_r_align_gocue = [];
n_mov_onset_trials = [];



for d = 1:length(reach_ds)
    
    
    
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
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
    ssc1 = stdata{1}.target{end}.neural_data.dim_red.st_scores(:,1:mani_dims,:);
    ssc2 = stdata{2}.target{end}.neural_data.dim_red.st_scores(:,1:mani_dims,:);
    
    emg1 = stdata{1}.target{end}.emg_data.emg;
    emg2 = stdata{2}.target{end}.emg_data.emg;
    
    
    
    % -----------------------------------------------------------------
    % Preallocate to store results
    mov_onsets = nan(2,size(ssc1,3));
    
    % nbr. of EMGs
    n_emgs = size(emg1,2);
    
    % Figure to plot each trial
    if plot_single_trial_flg
        stfig = figure('name','Single trial EMGs'); hold on;
    end
    
    
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    % Find movement onset for each trial
    for i = 1:size(ssc1,3)
        
        
        % --------------------------------------------------------------------------------
        % Compute SNR to choose subset of muscles that will
        % be used to estimate movement onset for this trial
        
        % SNR as peak2peak (because the EMGs are rectified
        % and normalized already)
        snr_emg1 = peak2peak(abs(emg1(:,:,i)),1);
        snr_emg2 = peak2peak(abs(emg2(:,:,i)),1);
        
%         % SNR as CoV
%         snr_emg1 = std(squeeze(emg1(:,:,i)),0,1)./mean(squeeze(emg1(:,:,i)),1);
%         snr_emg2 = std(squeeze(emg2(:,:,i)),0,1)./mean(squeeze(emg2(:,:,i)),1);
%         
%         % SNR as in Matt's paper (Perich & Miller EBR '17):
%         % peak2peak / 2*SD
%         snr_emg1 = peak2peak(emg1(:,:,i),1)./(2*std(emg1(:,:,i),0,1));
%         snr_emg2 = peak2peak(emg2(:,:,i),1)./(2*std(emg2(:,:,i),0,1));
        
        
%         % --------------------------------------------------------------------------------
%         % Plot EMGs and COV
%         
%         for l = 1:size(emg1,2)
%             lg1{l} = num2str(snr_emg1(l));
%             lg2{l} = num2str(snr_emg2(l));
%         end
%         
%         figure('name','EMGs task 1', 'units','normalized','outerposition',[.2 .2 .6 .6])
%         subplot(121); plot(squeeze(emg1(:,1:floor(n_emgs/2),i)),'linewidth',2)
%         legend(lg1(1:floor(n_emgs/2)),'Location','NorthEast'), legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
%         subplot(122); plot(squeeze(emg1(:,floor(n_emgs/2)+1:end,i)),'linewidth',2)
%         legend(lg1(floor(n_emgs/2)+1:end),'Location','NorthEast'), legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
%         
%         figure('name','EMGs task 2', 'units','normalized','outerposition',[.2 .2 .6 .6])
%         subplot(121); plot(squeeze(emg2(:,1:floor(n_emgs/2),i)),'linewidth',2)
%         legend(lg2(1:floor(n_emgs/2)),'Location','NorthEast'), legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
%         subplot(122); plot(squeeze(emg2(:,floor(n_emgs/2)+1:end,i)),'linewidth',2)
%         legend(lg2(floor(n_emgs/2)+1:end),'Location','NorthEast'), legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
        
        
        % --------------------------------------------------------------------------------
        % "Summed EMG": the movement-related signal that we
        % use to estimate movement onset
        
        emgs_use1 = snr_emg1 > snr_th_emg;
        emgs_use2 = snr_emg2 > snr_th_emg;
        
        nbr_emgs_use1 = sum(emgs_use1);
        nbr_emgs_use2 = sum(emgs_use2);
        
        % if there are at least 1 sufficiently modulated
        % EMG, find movement onset
        if nbr_emgs_use1 > 0 && nbr_emgs_use2 > 0
            
            if plot_single_trial_flg
                disp(['Trial: ' num2str(i)]);
                disp(['Nbr EMGs use 1: ' num2str(sum(emgs_use1))]);
                disp(['Nbr EMGs use 2: ' num2str(sum(emgs_use2))]);
            end
            
            % 0. "Summed EMG": sum squares EMGs with enough
            % SNR. EMGs may be different across tasks
            tot_emg1 = squeeze(sum(emg1(:,emgs_use1,i),2).^2);
            tot_emg2 = squeeze(sum(emg2(:,emgs_use2,i),2).^2);
            
            % 1. Postprocess the summed EMG
            % low pass filter and normalize to 0-1
            lp1 = filt_nodelay(tot_emg1,1/bin_size,fc_lp_emg);
            lp2 = filt_nodelay(tot_emg2,1/bin_size,fc_lp_emg);
            
            lp1 = lp1/max(lp1);
            lp2 = lp2/max(lp2);
            
            
            % --------------------------------------------------
            % If there's a lot of EMG at the beginning of
            % the trial exclude it
            if sum(lp1(1:20)>max_emg_beginning_trial) + sum(lp2(1:20)>max_emg_beginning_trial) == 0
                
                
                % find minima
                [min1,xmin1] = findpeaks(-lp1);
                [min2,xmin2] = findpeaks(-lp2);
                
                % find maxima
                [mx1, xmx1] = findpeaks(lp1);
                [mx2, xmx2] = findpeaks(lp2);
                
                % find maximum of maxima
                [supermx1, xsupermx1] = max(mx1);
                xsupermx1 = xmx1(xsupermx1);
                
                [supermx2, xsupermx2] = max(mx2);
                xsupermx2 = xmx2(xsupermx2);
            end
            
            % --------------------------------------------------
            % If one of the EMGs has no minima, exclude the trial
            if ~isempty(xmin1) && ~isempty(xmin2)
                
                
                % plot single trials?
                if plot_single_trial_flg
                    clf
                    figure(stfig), hold on
                    plot(lp1,'c');
                    for k = 1:length(min1); plot(xmin1(k),-min1(k),'.c','markersize',20); end
                    for k = 1:length(mx1); plot(xmx1(k),mx1(k),'oc','markersize',8); end
                    plot(lp2,'color',[.7 .7 .7]);
                    for k = 1:length(min2); plot(xmin2(k),-min2(k),'.','color',[.7 .7 .7],'markersize',20); end
                    for k = 1:length(mx2); plot(xmx2(k),mx2(k),'o','color',[.7 .7 .7],'markersize',8); end
                    title(['Trial: ' num2str(i)]);
                    set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
                    xlabel('Sample nbr.'), ylabel('"Summed EMG"')
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % FIND MOVEMENT ONSET
                
                % 2. find minima right before the maximum => this is
                % our movement onset
                % -- Note that there are a couple tricks to
                % prevent errors in weird trials
                minaftermax1 = find(xmin1>xsupermx1,1);
                minaftermax2 = find(xmin2>xsupermx2,1);
                
                % FIX #1: if there's no minima before the
                % absoluyte maximum, just take the one after it
                if minaftermax1 > 1
                    mo1 = xmin1(find(xmin1>xsupermx1,1)-1);
                else
                    mo1 = xmin1(find(xmin1>xsupermx1,1));
                end
                if minaftermax2 > 1
                    mo2 = xmin2(find(xmin2>xsupermx2,1)-1);
                else
                    mo2 = xmin2(find(xmin2>xsupermx2,1));
                end
                
                % FIX #2: minima don't come after maxima
                if isempty(mo1), mo1 = xmin1(end); end
                if isempty(mo2), mo2 = xmin2(end); end
                
                
                % get y coordinate of movement onsets
                ymo1 = lp1(mo1);
                ymo2 = lp2(mo2);
                
                % FIX #3: if there are two bursts of activity
                % and movement onset is found between bursts,
                % take the mimum before it
                if ymo1 > max_emg_mov_onset
                    if mo1 ~= xmin1(1)
                        mo1 = xmin1(find(xmin1==mo1)-1);
                        ymo1 = lp1(mo1);
                    end
                end
                if ymo2 > max_emg_mov_onset
                    if mo2 ~= xmin2(1)
                        mo2 = xmin2(find(xmin2==mo2)-1);
                        ymo2 = lp2(mo2);
                    end
                end
                
                % FIX #4: if there is a max that is quite large
                % before movement onset, find the minimum
                % before it
                emg_burst_before1 = find(mx1(xmx1 < mo1) > max_emg_mov_onset);
                while ~isempty( emg_burst_before1 ) && mo1 > xmin1(1)
                    % find position EMG burst before current
                    % movement onset estimate
                    if mo1 ~= xmin1(1)
                        mo1 = xmin1(find(xmin1==mo1)-1);
                    end
                    emg_burst_before1 = find(mx1(xmx1 < mo1) > max_emg_mov_onset);
                end
                emg_burst_before2 = find(mx2(xmx2 < mo2) > max_emg_mov_onset);
                while ~isempty( emg_burst_before2 ) && mo2 > xmin2(1)
                    if mo2 ~= xmin2(1)
                        mo2 = xmin2(find(xmin2==mo2)-1);
                    end
                    emg_burst_before2 = find(mx2(xmx2 < mo2) > max_emg_mov_onset);
                end
                
                % FIX #5: sometimes (because of the filtering)
                % the minimum is a negative EMG. This seems to
                % always be preceded by a positive peak. So
                % take the minimum before it
                if ymo1 < -0.05
                    prev_min1 = find(xmin1==mo1)-1;
                    if prev_min1 > 0
                        mo1 = xmin1(prev_min1);
                    end
                end
                if ymo2 < -0.05
                    prev_min2 = find(xmin2==mo2)-1;
                    if prev_min2 > 0
                        mo2 = xmin2(prev_min2);
                    end
                end
                
                
                % update y coordinate of movement onsets
                ymo1 = lp1(mo1);
                ymo2 = lp2(mo2);
                
                % un-invert for plotting
                min1 = -min1;
                min2 = -min2;
                
                
                % ------------------------------------------------------------------------------
                % Store results & Plot
                
                mov_onsets(:,i) = [mo1; mo2];
                
                if plot_single_trial_flg
                    figure(stfig), hold on
                    plot(mo1,ymo1,'c','marker','*','linewidth',3,'markersize',20)
                    plot(mo2,ymo2,'color',[.7 .7 .7],'marker','*','linewidth',3,'markersize',20)
                    pause(sec_pause);
                    clf
                end
            end
        else
            disp(['Excluding trial ' num2str(i)]);
        end
    end
    
    
    
    
    %% ------------------------------------------------------------------------
    % --------------------------------------------------------------------------
    % Exclude trials in which movement onset is too close to the end or the
    % beginning
    
    
    % to count stuff -- not necessary
    ctr_discard_beg = 0;
    ctr_discard_end = 0;
    ctr_discard_separated = 0;
    ctr_too_short = 0;
    
    for tr = 1:size(mov_onsets,2)
        
        % 1) too close to the beginning
        if mov_onsets(1,tr) < min_nbr_bins_beginning_mov_onset
            mov_onsets(:,tr) = nan(2,1);
            ctr_discard_beg = ctr_discard_beg + 1;
        elseif mov_onsets(2,tr) < min_nbr_bins_beginning_mov_onset
            mov_onsets(:,tr) = nan(2,1);
            ctr_discard_beg = ctr_discard_beg + 1;
        end
        
        % 2) too close to the end
        if ( size(emg1,1) - mov_onsets(1,tr) + 1 ) < min_nbr_bins_mov_onset_end
            mov_onsets(:,tr) = nan(2,1);
            ctr_discard_end = ctr_discard_end + 1;
        elseif ( size(emg2,1) - mov_onsets(2,tr) + 1 ) < min_nbr_bins_mov_onset_end
            mov_onsets(:,tr) = nan(2,1);
            ctr_discard_end = ctr_discard_end + 1;
        end
        
        % 3) too different across trials
        if abs( mov_onsets(1,tr) - mov_onsets(2,tr) ) > max_diff_mov_onset
            mov_onsets(:,tr) = nan(2,1);
            ctr_discard_separated = ctr_discard_separated + 1;
        end
        
        % 4) not including the number of bins we defined
        if ( size(emg1,1) - mov_onsets(1,tr) - n_bins_trial + trial_begins ) < 0
            mov_onsets(:,tr) = nan(2,1);
            ctr_too_short = ctr_too_short + 1;
        elseif ( size(emg2,1) - mov_onsets(2,tr) - n_bins_trial + trial_begins ) < 0
            mov_onsets(:,tr) = nan(2,1);
            ctr_too_short = ctr_too_short + 1;
        end
    end
    
    
    
    
    %% ------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Do CCA
    
    
    % -------------------------------------------------------------------------
    % Populate matrix with latent variables (concatenated trials)
    
    latvars1 = [];
    latvars2 = [];
    
    for tr = 1:size(mov_onsets,2)
        
        if ~isnan(mov_onsets(1,tr))
            
            trial1_beg = mov_onsets(1,tr) - trial_begins;
            %         if trial1_beg < 0, disp(['trial1_beg < 0 for tr = ' num2str(tr)]); end
            trial1_end = trial1_beg + n_bins_trial - 1;
            %         if trial1_end > size(ssc1,1), disp(['trial1_end > size(ssc1,1) for tr = ' num2str(tr)]); end
            latvars1 = [latvars1; ssc1(trial1_beg:trial1_end,:,tr)];
            
            trial2_beg = mov_onsets(2,tr) - trial_begins;
            trial2_end = trial2_beg + n_bins_trial - 1;
            latvars2 = [latvars2; ssc2(trial2_beg:trial2_end,:,tr)];
        end
    end
    
    
    [~,~,r] = canoncorr(latvars1,latvars2);
    
    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % compare to our old way of doing it...
    
    stdata_oldway = datasets{reach_ds(d)}.stdata;
    
    % trim original data window to same length as the current one?
    if trim_original_win_flg
        win_oldway = [proj_params.time_win(reach_ds(d),1), proj_params.time_win(reach_ds(d),1)+bin_size*n_bins_trial];
    else
        win_oldway = proj_params.time_win(reach_ds(d),:);
    end
    
    % 1) equalize trial duration across all tasks
    stdata_oldway = equalize_single_trial_dur( stdata_oldway, ...
        'time_win', win_oldway );
    
    % 2) equalize number of trials for all targets of a given task
    for i = 1:length(stdata_oldway)
        stdata_oldway{i} = equalize_nbr_trials_p_target( stdata_oldway{i} );
    end
    
    % 3) equalize number of trials across tasks
    stdata_oldway = equalize_nbr_trials_across_tasks( stdata_oldway, 'all_conc' );
    
    latvars1_oldway = stdata_oldway{1}.target{end}.neural_data.dim_red.scores(:,1:mani_dims);
    latvars2_oldway = stdata_oldway{2}.target{end}.neural_data.dim_red.scores(:,1:mani_dims);
    
    
    [~,~,r_oldway] = canoncorr(latvars1_oldway,latvars2_oldway);
    
    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Store values
    all_r = [all_r; r];
    all_r_align_gocue = [all_r_align_gocue; r_oldway];
    
    n_mov_onset_trials = [n_mov_onset_trials, sum(~isnan(mov_onsets(1,:)))]; 
    
    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Plot for this session
    
    if plot_each_session_flg
        res_ds = figure;
        plot(r,'k','linewidth',2)
        hold on,plot(r_oldway,'linewidth',2,'color',[.6 .6 .6])
        set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
        xlabel('Neural mode'), ylabel('Cross-task canonical correlation')
        legend('align to mov. onset','align to go signal'), legend boxoff
        title([datasets{reach_ds(d)}.monkey ' ' datasets{reach_ds(d)}.date(1:end-4)]);
        ylim([0 1]); title(['N. trials: ' num2str(n_mov_onset_trials(d))])
        pause;
        close(res_ds);
    end
    
    
    %% ------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Do new shuffling control: shuffle neural unit firing across channels
    % and trials
    
    if do_control_flg
        
    end
end





%% ------------------------------------------------------------------------
%% SUMMARY PLOTS


colors = parula(numel(reach_ds)+1);



% -------------------------------------------------------------------------
% 1. Plot CC aligning to movement onset vs CC aligning to go cue

figure,hold on
for d = 1:length(reach_ds)
    plot(all_r(d,:),'linewidth',2,'color',colors(d,:));
    plot(all_r_align_gocue(d,:),'-.','linewidth',2,'color',colors(d,:));
end
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlabel('Neural mode'), ylabel('Cross-task canonical correlation')
legend('align to mov. onset','align to go cue'), legend boxoff



% -------------------------------------------------------------------------
% 2. Plot ratio CCs between methods

lg_pooled = {};
for d = 1:length(reach_ds)
%     lg_pooled{d} = [datasets{reach_ds(d)}.monkey(1) ' - ' datasets{reach_ds(d)}.date(1:end-4) ' - n=' num2str(n_mov_onset_trials(d))];
    lg_pooled{d} = [datasets{reach_ds(d)}.monkey(1) ' - ' datasets{reach_ds(d)}.date(1:end-4)];
end

figure,hold on
for d = 1:length(reach_ds)
    plot(all_r(d,:)./all_r_align_gocue(d,:),'linewidth',2,'color',colors(d,:));
    %plot(all_r_align_gocue(d,:),'-.','linewidth',2,'color',colors(d,:));
end
plot([1 mani_dims],[1 1],'-.','color',[.7 .7 .7],'linewidth',3)
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
xlabel('Neural mode'), ylabel('CC aligning mov. onset / CC aligning go cue')
xlim([0 mani_dims+5])
title(['Nbr bins: ' num2str(n_bins_trial) ' - trial begins (wrt mov onset): ' num2str(-trial_begins)]);
legend(lg_pooled), legend boxoff
