%
% Compute SNR of the EMGs and neural data --to respond to Reviewer 3
% Comment 1b 
%
% The computation of the SNR is from (Perich & Miller, 2017)

clc;

ds_to_use = [1:3 7:9];

% define datasets
wrist_ds = [1:3 7:9];
reach_ds = [4:6 10:11];


% define manifold dimensionality, to know how many neural modes to keep
mani_dim = 12;



%% ------------------------------------------------------------------------
% Arrange the data for the analysis


% variable names
vars = {'units','emgs','modes'};

% Define matrices to store data
for v = 1:length(vars)
    
    p2p.(vars{v}) = [];
    mn.(vars{v}) = [];
    sd_inter_trial.(vars{v}) = [];
    wrist_datum.(vars{v}) = []; % flag to say if it's wrist or reach to grasp
end
 


% Fill them with the data goods!
for d = 1:length(ds_to_use)
   
    % dataset ptr
    t_ds = ds_to_use(d);
    
    for t = 1:length(datasets{t_ds}.labels)
        
        % ------------------------------------------------------------
        % Get Mean and Peak to Peak value of the trial-averaged data
        
        % Units
        mn.units = [mn.units; mean(datasets{t_ds}.stdata{t}.target{end}.neural_data.smoothed_fr_mn,1)'];
        p2p.units = [p2p.units; peak2peak(datasets{t_ds}.stdata{t}.target{end}.neural_data.smoothed_fr_mn,1)'];
        
        % EMGs
        mn.emgs = [mn.emgs; mean(datasets{t_ds}.stdata{t}.target{end}.emg_data.mn,1)'];
        p2p.emgs = [p2p.emgs; peak2peak(datasets{t_ds}.stdata{t}.target{end}.emg_data.mn,1)'];
        
        % neural modes
        mn.modes = [mn.modes; mean(datasets{t_ds}.stdata{t}.target{end}.neural_data.dim_red.st_scores_mn(:,1:mani_dim),1)'];
        p2p.modes = [p2p.modes; peak2peak(datasets{t_ds}.stdata{t}.target{end}.neural_data.dim_red.st_scores_mn(:,1:mani_dim),1)'];
        
        
        % ------------------------------------------------------------
        % Calculate SD of the mean-subtracted data
        
        % pre-allocate matrices to store the inter-trial fluctuations
        t_units_inter_trial = [];
        t_emgs_inter_trial = [];
        t_modes_inter_trial = [];
        
        for tg = 1:length(datasets{t_ds}.stdata{t}.target)-1
           
            % -----------------
            % Units
            % repeat the trial-average activity
            t_units_mn = repmat(datasets{t_ds}.stdata{t}.target{tg}.neural_data.smoothed_fr_mn,...
                            size(datasets{t_ds}.stdata{t}.target{tg}.neural_data.smoothed_fr,3),1);
            
            % and subtract from the single trial data to get the
            % inter-trial fluctuations
            t_units_inter_trial = [t_units_inter_trial; datasets{t_ds}.stdata{t}.target{tg}.neural_data.conc_smoothed_fr - t_units_mn];
            
            % -----------------
            % Neural modes
            t_modes_mn = repmat(datasets{t_ds}.stdata{t}.target{tg}.neural_data.dim_red.st_scores_mn(:,1:mani_dim),...
                            size(datasets{t_ds}.stdata{t}.target{tg}.neural_data.dim_red.st_scores,3),1);
            
            t_modes_inter_trial = [t_modes_inter_trial; datasets{t_ds}.stdata{t}.target{tg}.neural_data.dim_red.scores(:,1:mani_dim) - t_modes_mn];
            
            % -----------------
            % EMGs
            t_emgs_mn = repmat(datasets{t_ds}.stdata{t}.target{tg}.emg_data.mn,...
                            size(datasets{t_ds}.stdata{t}.target{tg}.emg_data.emg,3),1);
                        
            t_emgs_inter_trial = [t_emgs_inter_trial; datasets{t_ds}.stdata{t}.target{tg}.emg_data.conc_emg - t_emgs_mn];
        end
        
        % -----------------
        % Compute the SD of the inter-trial fluctuations
        
        sd_inter_trial.units = [sd_inter_trial.units; std(t_units_inter_trial,0,1)'];
        sd_inter_trial.emgs = [sd_inter_trial.emgs; std(t_emgs_inter_trial,0,1)'];
        sd_inter_trial.modes = [sd_inter_trial.modes; std(t_modes_inter_trial,0,1)'];
        
        
        % ------------------------------------------------------------
        % Wrist data? -- flag that says if it's reach or wrist
        
        n.units = length(datasets{t_ds}.neural_chs);
        n.modes = mani_dim;
        n.emgs = length(datasets{t_ds}.chosen_emgs);
        
        if ismember(ds_to_use(d),wrist_ds)
            wrist_flag = 1;
        else
            wrist_flag = 0;
        end
        
        for v = 1:length(vars)
            wrist_datum.(vars{v}) = [wrist_datum.(vars{v}); wrist_flag*ones(n.(vars{v}),1)];
        end
    end
end



%% ------------------------------------------------------------------------
% Compute SNR
%
% SNR (Perich & Miller, 2017): Waveform peak-to-peak value divided by two
% times the SD of the mean-subtracted data. If the distributions are very
% different add noise to the neural modes


for v = 1:length(vars)
    
    snr.(vars{v}) = p2p.(vars{v})./(2*sd_inter_trial.(vars{v}));
end


% Do a histogram and compute summary stats

hist_bin = .2;
x_hist = 0:hist_bin:10+hist_bin;

for v = 1:length(vars)
   
    hist_snr.(vars{v}) = histcounts(snr.(vars{v}),x_hist)/numel(snr.(vars{v}));
    mn_snr.(vars{v}) = mean(snr.(vars{v}));
    sd_snr.(vars{v}) = std(snr.(vars{v}));
end


% Do paired t-tests to see if the distributions are different
[~, p_emgs_modes] = ttest2(snr.emgs,snr.modes);
[~, p_emgs_units] = ttest2(snr.emgs,snr.units);
disp(['Probably distribs. EMGs and Neural Modes SNR is similar (t-test): ' num2str(p_emgs_modes)]);
disp(['Probably distribs. EMGs and Units SNR is similar (t-test): ' num2str(p_emgs_units)]);



%% ------------------------------------------------------------------------
% PLOT

col.units = [0 204 102]/255;
col.emgs = [247 152 80]/255;
col.modes = [157 157 182]/255;

figure,hold on
bar(x_hist(1:end-1),hist_snr.units*100,'FaceColor',col.units)
bar(x_hist(1:end-1),hist_snr.emgs*100,'FaceColor',col.emgs), alpha(gca,.5)
bar(x_hist(1:end-1),hist_snr.modes*100,'FaceColor',col.modes), alpha(gca,.5)
% plot errorbars
yl = ylim;
for v = 1:length(vars)
    plot(mn_snr.(vars{v}),yl(2)+(v-1.5),'.','markersize',20,'color',col.(vars{v}))
    plot([mn_snr.(vars{v})-sd_snr.(vars{v}),mn_snr.(vars{v})+sd_snr.(vars{v})],[yl(2)+(v-1.5) yl(2)+(v-1.5)],'linewidth',1.5,'color',col.(vars{v}))
end
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Signal-to-noise ratio'), ylabel('Percentage (%)')
legend(['Units, n=' num2str(numel(mn.units))],['EMGs, n=' num2str(numel(mn.emgs))],['Modes, n=' num2str(numel(mn.modes))]),legend boxoff
