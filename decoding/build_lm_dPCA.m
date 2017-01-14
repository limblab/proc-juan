%
% Build linear model to predict EMG/force/kinematics using dPCs as inputs.
% Note: This function uses trial-averaged data !!!
%
%
%   dPCA_lm = build_lm_dPCA( dPCA_data, st_data, varargin )
%
%
% Inputs (opt)          : [default]
%   dPCA_data           : dPCA_results data struct, generated with
%                           call_dPCA.m
%   st_data             : single_trial_data struct
%   (neuron_emg_delay)  : [0.04] neuron to EMG lag (>0 => neurons lead)
%   (margs)             : [2] marginalizations (dPCs that depend on certain
%                           covariates) used as inputs. With the current
%                           version of call_dPCA they are: 1) 'task'; 2)
%                           'target'; 3) 'time'; 4) 'task/target interact'
% 

function dPCA_lm = build_lm_dPCA( dPCA_data, st_data, varargin )


if nargin >= 3
    neuron_emg_delay    = varargin{1};
else
    neuron_emg_delay    = 0.04;
end
if nargin >= 4
    margs               = varargin{2};
else
    margs               = 2; 
end
if nargin >= 5
    xval_yn             = varargin{3};
else
    xval_yn             = true;
end
if nargin == 6
    xval_plot_yn        = varargin{4};
else
    xval_plot_yn        = false;
end


delay_bins      = neuron_emg_delay / st_data{1}.target{1}.bin_size;


% get the dPCA latent variables
lvars           = dPCA_data.lat_vars_mn;


% -------------------------------------------------------------------------
% 1. Some preliminary stuff for the EMG data -- copied from call_dPCA

% make the target averaged responses for each task have equal length.
% This is not done in single_trial_analysis.m, where single trial
% duration is only equalized for each task
st_data        	= equalize_single_trial_dur( st_data );

% get rid of the last target, which is all the concatenated targets
for i = 1:length(st_data)
    st_data{i}.target(end) = [];
end


% -------------------------------------------------------------------------
% 2. Preprare EMG data

% create a 4D matrix with the EMG/forces/kinematics organized as in the
% dPCA matrix: 
% M is the number of muscles/forces/kinematic signals
% K is the number of tasks
% G is the number of targets
% T is the number of time points
M               = numel(st_data{1}.target{1}.emg_data.emg_names);
K               = length(st_data);
nbr_tgts_p_task = cell2mat( cellfun(@(x) length(x.target), st_data,'UniformOutput',false) );
G               = min(nbr_tgts_p_task);
T               = size(st_data{1}.target{1}.emg_data.mn,1);

% a) preallocate and fill the matrix of average EMGs
emg_mn          = nan(M,K,G,T);
% ·········································································
% In the 1D tasks, targets 1 to 6 go from left to right, but in the 2D
% tasks targets are ordered as follows: 5, 7, 8, 6, 4, 2, 1, 3 --beginning
% at 12 o'clock and going clockwise. They will be paired as 1D/2D: 1/1,
% 2/2, 3/3, 4/6, 5/7, 6, 8 (as defined in target_order)
% ·········································································
target_order    = [1 2 3 4 5 6; 1 2 3 6 7 8]; % row 1: 1D task; row 2: 2D task

for k = 1:K
    for g = 1:G
        if length(st_data{k}.target) == G
            tgt = target_order(1,g);
        elseif length(st_data{k}.target) >= G
            tgt = target_order(2,g);
        end
        emg_mn(:,k,g,:) = st_data{k}.target{tgt}.emg_data.mn';
    end    
end


% b) take the neuron to EMG delay into account

if delay_bins >= 0
    emg_mn      = emg_mn(:,:,:,1+delay_bins:T);
elseif delay_bins < 0
    emg_mn      = emg_mn(:,:,:,1:T-delay_bins);
end

T_new           = size(emg_mn,4);


% c) create a matrix that has the EMG for all the targets of each task
% ordered sequentially
emg_lm          = zeros(M,K*G*T_new);
for m = 1:M
    for k = 1:K
        for g = 1:G
            emg_lm( m, (1:T_new) + (g-1)*T_new + (k-1)*G*T_new ) = ...
                squeeze(emg_mn(m,k,g,:))';
        end
    end
end


% -------------------------------------------------------------------------
% 2. Preprare dPCA (neural) data

lv              = dPCA_data.lat_vars;
N               = size(lv,1);

% a) take the neuron to EMG delay into account

if delay_bins > 0
    lv          = lv(:,:,:,1:T-delay_bins);
else
    lv          = lv(:,:,:,1+delay_bins:T);
end

% b) create a matrix that has the latent variables for all the targets of
% each task ordered sequentially
lv_lm           = zeros(N,K*G*T_new);
for n = 1:N
    for k = 1:K
        for g = 1:G
            lv_lm( n, (1:T_new) + (g-1)*T_new + (k-1)*G*T_new ) = ...
                squeeze(lv(n,k,g,:))';
        end
    end
end


% -------------------------------------------------------------------------
% 3. Build linear model

y               = emg_lm';
X               = lv_lm';

% choose the latent variables that we want. By default, the
% marginalizations are 'task', 'target', 'time', 'task/target interact'
lvs_marg        = [];
for i = 1:length(margs)
    lvs_marg    = [lvs_marg , find(dPCA_data.which_marg==margs(i))];
end
X               = X(:,lvs_marg);


if xval_yn
    [W, Wos, stats] = build_lm( X, y, true, length(X)/K, xval_plot_yn );
else
    [W, Wos, stats] = build_lm( X, y, false );
end


% -------------------------------------------------------------------------
% 4. Return vars

dPCA_lm.W       = W;         
dPCA_lm.stats   = stats;
dPCA_lm.margs   = margs;
dPCA_lm.data.emg    = y;
dPCA_lm.data.pred   = W*X'+repmat(Wos,1,size(X,1));
