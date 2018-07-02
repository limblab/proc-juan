%
% Compare predictions of different EMG marginaliaztions from target-related
% neural modes
%


% SAME PARAMETERS AS IN COMP_DECODERS_DPCA

% crossvalidate the predicitons?
xval_yn     = true;

% neural manifold dimensionality
mani_dim    = 12;

% emg manifold dimensionality
emg_mani_dim = 4;

% number of dPCs per marginalization to use
dpcs_marg   = 2;

% use s.e.m. rather than s.d. for plotting
use_sem_yn  = false;

% fold size (s)
fold_size   = 30;

% order polynomial static non-linearity
poly_order  = 2;

% marginalization legends for plotting fits -- THIS DEFINES DATA IN ALL THE FIGURES
marg_leg    = {'task','target','time','task/tgt int'};

% combinations of marginalizations that will be compared
% marg_set    = {1, 2, 3, 4, [2 4], [1 3], [2 3], [1 2 4], [1 2 3], [1 3 4], [1 2 3 4]};
marg_set    = {1, 2, 3, 4, [1 2 3 4]};

% datasets that will be analyzed
D           = [1:3 7:9]; % Jaco and Jango 1D / 2D datasets

% use only good EMGs? --- THIS NEEDS TO BE BUILT IN BUILD_DP_DECODER, in
% the mean time is used for plotting and for the stats
good_emgs_only = true;



for d = 1:length(D)
    if good_emgs_only
        emg_chs{d} = datasets{D(d)}.chosen_emgs;
    else
        emg_chs{d} = 1:length(datasets{D(d)}.binned_data{1}.emgguide);
    end
end


num_margs   = length(marg_set);



%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 0. Do dPCA of the neural data, if the data are not avaialble

if exist('dPCA_results','var')
    dPCA_data = dPCA_results;
end
if ~exist('dPCA_data','var')
    for d = 1:length(D)
        dPCA_data{d} = call_dPCA( datasets{D(d)}.stdata, mani_dim, false );
    end
end


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 1. Do dPCA of the EMGs

% for d = 1:5%length(D)
%     dPCA_emg{d} = call_dPCA_emgs( datasets{D(d)}.stdata, emg_mani_dim, true);
% end


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 2. Predict EMGs for each marginalization
    
for d = 1:length(D)
    
    disp(['building decoders with chosen marginalizations for dataset ', num2str(D(d))]);
    disp('...');
    for m = 1:length(marg_set)-1
        %if m == 2, plot_flag = true; else plot_flag = false; end
        dPCA_fit{d,m} = build_dPC_decoder_new( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'margs',            marg_set{2}, ...
                            'nbr_dpcs_p_marg',  dpcs_marg, ... 
                            'xval_yn',          xval_yn, ...
                            'fold_size',        fold_size, ...
                            'smooth_spikes',    false, ...
                            'dec_offset',       true, ...
                            'poly_order',       poly_order, ...
                            'return_preds',     true, ...
                            'dec_p_task',       false, ...
                            'tgts_to_excl',     {[-1 -5 1 -7],[-1 7 1 5]}, ...
                            'plot_yn',          false, ...
                            'output',           'dpc_emg', ...
                            'dims_dpca_emg',    emg_mani_dim, ...
                            'marg_dpca_emg',    marg_set{m});
    end
end


%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 3. Predict EMGs with all dPCAs
    
for d = 1:length(D)
    
    disp(['building decoders with all dPCs marginalizations for dataset ', num2str(D(d))]);
    disp('...');
    %if m == 2, plot_flag = true; else plot_flag = false; end
    dPCA_fit_all_dPCs{d} = build_dPC_decoder_new( datasets{D(d)}, dPCA_data{d}, ...
                            'lags',             10, ...
                            'dpcs',             1:mani_dim, ...
                            'xval_yn',          xval_yn, ...
                            'fold_size',        fold_size, ...
                            'smooth_spikes',    false, ...
                            'dec_offset',       true, ...
                            'poly_order',       poly_order, ...
                            'return_preds',     true, ...
                            'dec_p_task',       false, ...
                            'tgts_to_excl',     {[-1 -5 1 -7],[-1 7 1 5]}, ...
                            'plot_yn',          false, ...
                            'output',           'dpc_emg', ...
                            'dims_dpca_emg',    emg_mani_dim, ...
                            'marg_dpca_emg',    2); 
end



%% -----------------------------------------------------------------------
% ------------------------------------------------------------------------
% 4. Pool results per marginalization across sessions (both monkeys
% together)
    

% R2 all muscles selected muscles for each marginalization
R2_ds_muscle_marg       = cell(length(D),length(marg_leg));

for d = 1:length(D)
    for m = 1:length(marg_leg)
        if ~isempty( dPCA_fit{d,m}.mfxval{1}.R2 )
            R2_ds_muscle_marg{d,m} = dPCA_fit{d,m}.mfxval{1}.R2(1,:);
        else
            R2_ds_muscle_marg{d,m} = nan(1,size(dPCA_fit{d,m}.mfxval{1}.R2,2));
        end
    end
end

% R2 EMGs for all selected muscles using all dPCs as predictors
R2_ds_muscle_marg_all_dPCs   = cell(length(D),1);

for d = 1:length(D)
    R2_ds_muscle_marg_all_dPCs{d} = dPCA_fit_all_dPCs{d}.mfxval{1}.R2(1,:);
end


% POOL R2s ONTO MATRICES


for m = 1:length(marg_leg)
    R2_per_marg_mtrx(:,m) = [R2_ds_muscle_marg{:,m}];
end

R2_all_dPCs_mtrx = cell2mat(R2_ds_muscle_marg_all_dPCs')';

norm_R2_marg = R2_per_marg_mtrx ./ R2_all_dPCs_mtrx;

