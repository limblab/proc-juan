%
%  
%

clear all, %close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters / definitions
monkey          = 'mihili'; % 'mihili'; 'chewie'

spiking_inputs  = {'M1_spikes'};

% The behavioral variable to be predicted
dec_var         = 'vel';
% Decoder inputs
dec_input       = 'aligned_data'; % 'aligned_data'; 'raw_data'

% When we start looking at each trial
idx_start       = {'idx_movement_on',0};
% When we stop looking at each trial --if empty, the duration of the shortest trial % {'idx_go_cue',60}
idx_end         = {''};
% Neural modes that will be aligned
mani_dims       = 1:4;
% Downsampling rate (nbr of bins that will be combined)
n_bins          = 5;
% N bins for decoder history
hist_bins       = 2;


files_chewie = { ...
                'Chewie_CO_VR_2016-09-12.mat', ...
                'Chewie_CO_VR_2016-09-14.mat', ...
                'Chewie_CO_FF_2016-09-15.mat', ...
                'Chewie_CO_FF_2016-09-19.mat', ...
                'Chewie_CO_FF_2016-09-21.mat', ...
                'Chewie_CO_FF_2016-10-05.mat', ...
                'Chewie_CO_VR_2016-10-06.mat', ...
                'Chewie_CO_FF_2016-10-07.mat', ...
                'Chewie_CO_FF_2016-10-11.mat', ...
                'Chewie_CO_FF_2016-10-13.mat' ...
                };

files_mihili = { ...
                'Mihili_CO_FF_2014-02-03.mat', ...
                'Mihili_CO_FF_2014-02-17.mat', ...
                'Mihili_CO_FF_2014-02-18.mat', ...
                'Mihili_CO_VR_2014-03-03.mat', ...
                'Mihili_CO_VR_2014-03-04.mat', ...
                'Mihili_CO_VR_2014-03-06.mat', ...
                'Mihili_CO_FF_2014-03-07.mat' ...
                };
 

% override monkey name and go to path
here        = pwd;
switch monkey
    case 'chewie'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/Chewie')
        files = files_chewie;
    case 'mihili'
        cd('/Users/juangallego/Documents/NeuroPlast/Data/Mihili')
        files = files_mihili;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data, unsorting the channels. Baseline trials only

master_td   = loadTDfiles(  files, ...
                            {@mergeSortedNeurons}, ...
                            {@getTDidx,'epoch','BL'}, ...
                            {@getTDidx,'result','R'}, ...
                            {@sqrtTransform,spiking_inputs}, ...
                            {@smoothSignals,struct('signals',spiking_inputs,'calc_fr',true,'kernel_SD',0.05)}, ...
                            {@getPCA,struct('signals',spiking_inputs)}, ...
                            {@binTD,n_bins}...
                            );

                        
% go back to where you were path-wise
cd(here);


% get some meta info about the data

% get the sessions
sessions    = unique({master_td.date});

% get the targets
targets     = unique(cell2mat({master_td.target_direction}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim the trials

% get the minimum trial duration, if needed for idx_end
if isempty(idx_end{1})
    min_mov_duration = min( arrayfun( @(x) x.idx_trial_end-x.idx_movement_on, master_td) );
    idx_end = {'idx_movement_on',min_mov_duration};
end    

master_td   = trimTD( master_td, idx_start, idx_end );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep the same number of trials for all the targets and sessions --this is
% important to align the neural mode dynamics with CCA

master_td   = equalNbrTrialsSessions(master_td);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do CCA of the neural trajectories across sessions, to align the dynamics
% and see how similar they are

% get all pairwise comparisons of sessions
comb_sessions   = nchoosek(1:length(sessions),2);

% do
disp('Aligning neural mode dynamics across days')
for c = 1:size(comb_sessions,1)
    
    % get the data
    trials1 = getTDidx( master_td, 'date', sessions{comb_sessions(c,1)} );
    trials2 = getTDidx( master_td, 'date', sessions{comb_sessions(c,2)} );
    
    %trials2 = trials2(randperm(length(trials2)));
                
    % compare dynamics with CCA
    cca_info(c) = compDynamics( master_td, 'M1_pca', trials1, trials2, mani_dims );
    
    % compare dynamics with good old forrelations
    corr_info(c) = corrDynamics( master_td, 'M1_pca', trials1, trials2, mani_dims );
end


% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build a decoder on one day, and test it in another. Do this for all pairs
% of sessions

disp('Testing decoder stability');
disp(['Decoding: ' dec_var]);

for c = 1:size(comb_sessions,1)

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add aligned neural mode dynamics to trial_data
    
    % get two TD structs, one per session that we are comparing. The code
    % will build a model on TD1 and test it on TD2
    [trials1, td1] = getTDidx(master_td,'date',sessions{comb_sessions(c,1)});
    [trials2, td2] = getTDidx(master_td,'date',sessions{comb_sessions(c,2)});
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SANTIY CHECK: the CC between the aligned & unaligned data is all 1 as
    % it should
    %
    x = get_vars(td1,{'M1_pca',mani_dims});
    y = get_vars(td2,{'M1_pca',mani_dims});
    u = cca_info(c).U;
    v = cca_info(c).V;
    [~,~,cc_xu] = canoncorr(x,u);
    [~,~,cc_yv] = canoncorr(y,v);
    if min(cc_xu) < .99, error(['CC between raw and aligned data ~1 for c = ' num2str(c)]); end
    if min(cc_yv) < .99, error(['CC between raw and aligned data ~1 for c = ' num2str(c)]); end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Separate concatenated neural mode dynamics into trials
    bins_p_trial = size(td1(1).pos,1);
%     U = zeros(bins_p_trial,length(mani_dims),length(trials1));
%     V = zeros(bins_p_trial,length(mani_dims),length(trials1));
    for t = 1:length(trials1)
        be = bins_p_trial*(t-1)+1;
        en = bins_p_trial*t;
        td1(t).align_M1_pca = cca_info(c).U(be:en,:);
        td2(t).align_M1_pca = cca_info(c).V(be:en,:);
%         U(:,:,t) = cca_info(c).U(be:en,:);
%         V(:,:,t) = cca_info(c).V(be:en,:);
    end

%     % this adds a new field with the aligned data
%     for t = 1:length(trials1)
%         td1(t).align_M1_pca = U(:,:,t);
%         td2(t).align_M1_pca = V(:,:,t);
%     end
    
    % duplicate and shift --to have history in the model
    td1 = dupeAndShift(td1,{'align_M1_pca',hist_bins});
    td2 = dupeAndShift(td2,{'align_M1_pca',hist_bins});

    % build linear model
    mod_params.model_type = 'linmodel';
    mod_params.out_signals = dec_var;
    switch dec_input
        case 'aligned_data'
            mod_params.in_signals = {'align_M1_pca_shift',1:length(mani_dims)*hist_bins};
        case 'raw_data'
            td1 = dupeAndShift(td1,{'M1_pca',hist_bins});
            td2 = dupeAndShift(td2,{'M1_pca',hist_bins});
            mod_params.in_signals = {'M1_pca_shift',1:length(mani_dims)*hist_bins};
        otherwise
            error('wrong decoder input');
    end
    [td1, mod_info] = getModel(td1,mod_params);
    
    % compute R2 of within-day predictions
    Y = get_vars(td1,{dec_var,1:size(td1(1).(dec_var),2)});
    Yhat = get_vars(td1,{'linmodel_default',1:size(td1(1).(dec_var),2)});
      
    R2 = compute_r2(Y,Yhat,'corr');
    withinR2(c,:) = R2;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now test the model on the data from day 2
        
    [R2diff td2] = testModel(td2,mod_info);
    acrossR2(c,:) = R2diff;
end


figure,plot(unique(round(withinR2*100)/100,'rows'),'-o','linewidth',2),ylim([0 1])
set(gca,'TickDir','out','FontSize',14), box off
legend('X','Y','Location','SouthEast'),legend('boxoff')
ylabel('Model R^2'),xlabel('Session')
set(gca,'XTick',1:length(sessions)-1),xlim([0 length(sessions)])
switch dec_input
    case 'aligned_data'
        title(['Within-day predictions from an aligned manifold with D = ' num2str(length(mani_dims))]);
    case 'raw_data'
        title(['Within-day predictions from a manifold with D = ' num2str(length(mani_dims))]);
end


figure,plot(diff_days,acrossR2,'o','linewidth',2),ylim([0 1])
set(gca,'TickDir','out','FontSize',14), box off
legend('X','Y','Location','SouthEast'),legend('boxoff')
ylabel('Model stability over sessions R^2'),xlabel('Days between sessions')
