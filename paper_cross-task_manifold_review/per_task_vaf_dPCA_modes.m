%
% Compute per-task variance of the dPCA modes, as requested by Reviewer 2
% of the Nature Comms paper
% 


% choose manifold dimensionality
manifold_dim = 12;

% Use all wrist datasets for dPCA
dPCA_datasets = [1 2 3 7 8 9];

% What dPCA datasets correspond to each monkey
jaco_datasets = 1:3;
jango_datasets = 4:6;
% After pooling together all the results, define task indexes for each
% monkey
jaco_tasks = 1:9;
jango_tasks = 10:20;



% Do dPCA if the data aren't available in the WS
if ~exist('dPCA_results','var')
    
    for i = 1:length(dPCA_datasets)
        dPCA_results{i} = call_dPCA( datasets{dPCA_datasets(i)}.stdata, manifold_dim, false );
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we need to get the neural activity as we have projected it onto
% the manifold -- this code is copied directly from call_dPCA.m


for ds = 1:length(dPCA_datasets) % with respect to dPCA_datasets (!!!)

    
    % data for this session
    stdata = datasets{dPCA_datasets(ds)}.stdata;

    % neural units to use
    neural_chs = stdata{1}.target{1}.neural_data.neural_chs;

    % check if the dimensions in single_trial_data are consistent
    if size(stdata{1}.target{1}.neural_data.smoothed_fr,2) ~= numel(neural_chs)
        error('single trial data does not include all the neural channels')
    end

    % make the target averaged responses for each task have equal length.
    % This is not done in single_trial_analysis.m, where single trial
    % duration is only equalized for each task
    stdata = equalize_single_trial_dur( stdata );

    % get rid of the last target, which is all the concatenated targets
    for i = 1:length(stdata)
        stdata{i}.target(end) = [];
    end


    % ------------------------------------------------------------------------
    % 1. arrange the data (as described in dpca_demo)

    % N is the number of neurons
    % S is the number of conditions --> tasks in our case
    % D is the number of decisions --> targets in our case
    %       ToDo: so far we are choosing the min, but see if they can be different for each task
    % T is the number of time points --each trial should have the same duration
    % in time !!!
    N                   = numel(neural_chs);
    S                   = numel(stdata);
    D                   = min(cellfun(@(x) length(x.target), stdata));
    T                   = size(stdata{1}.target{1}.neural_data.fr,1);
    % max number of repetitions
    max_trial_num       = 0;
    for i = 1:S
        if max(cellfun(@(x) size(x.neural_data.fr,3), stdata{1}.target )) > max_trial_num
            max_trial_num = max(cellfun(@(x) size(x.neural_data.fr,3), stdata{1}.target ));
        end
    end


    % trial_num: N x S x D
    trial_num           = zeros(N,S,D);


    % firing_rates: N x S x D x T x max_trial_num -- populated with our
    % single_trial_data
    firing_rates        = nan(N,S,D,T,max_trial_num);
    % ·········································································
    % In the 1D tasks, targets 1 to 6 go from left to right, but in the 2D
    % tasks targets are ordered as follows: 5, 7, 8, 6, 4, 2, 1, 3 --beginning
    % at 12 o'clock and going clockwise. They will be paired as 1D/2D: 1/1,
    % 2/2, 3/3, 4/6, 5/7, 6, 8 (as defined in target_order)
    % ·········································································
    target_order        = [1 2 3 4 5 6; 1 2 3 6 7 8]; % row 1: 1D task; row 2: 2D task

    for n = 1:N
        for s = 1:S
            for d = 1:D
                if length(stdata{s}.target) == D
                    tgt = target_order(1,d);
                elseif length(stdata{s}.target) >= D
                    tgt = target_order(2,d);
                end
                trials_this = size(stdata{s}.target{tgt}.neural_data.smoothed_fr(:,n,:),3);
                firing_rates(n,s,d,:,1:trials_this) = ...
                    squeeze(stdata{s}.target{tgt}.neural_data.smoothed_fr(:,n,:));
            end
        end
    end


    % firing_rates_average: N x S x D x T -- these are PSTHs
    firing_rates_avg    = nanmean(firing_rates, 5);

    
    % ------------------------------------------------------------------------
    % Compute the variance explained by each dPCA mode by projecting the
    % neural data for each task onto the dimensions of the common manifold

    % get the basis of the dPCA manifold
    W = dPCA_results{ds}.W;
    V = dPCA_results{ds}.V;
    
    % Orthonormalize the basis
    Worth = orth(W);

    per_task_var = zeros(S,manifold_dim);
    per_task_var_per_marg = zeros(S,length(dPCA_results{ds}.marg_names));

    % do for each task
    for t = 1:S

        % Get trial-averaged FRs for this task
        fr_this = squeeze(firing_rates_avg(:,t,:,:));

        % Stitch together the average activity patterns for each target
        fr_this = permute(fr_this,[1 3 2]);
        conc_fr = reshape(fr_this, [size(fr_this,1), size(fr_this,2)*size(fr_this,3)]);

%         % Project the data onto the orthonormalized manifold dimenions
%         dPCA_modes_this = Worth'*conc_fr;

        
        % -----------------------------------------------------------------
        % Call Machens' dPCA explained variance function

        explVar = dpca_explainedVariance(conc_fr, W, V );
        
        % Compute how much of the total variance the dPCA modes explain for
        % each task -- DOES THIS MAKE SENSE? 
        per_task_var(t,:) = explVar.cumulativeDPCA;

        
        
%         % Add up the per task covariance per covariate (time, target, ...)
%         for c = 1:length(dPCA_results{ds}.marg_names)
%             per_task_var_per_marg(t,c) = sum( per_task_var(t,dPCA_results{ds}.which_marg==c) );
%         end
    end

    % ---------------------------------------------------------------------
    % Store the results
    
    dPCA_results{ds}.per_task_var = per_task_var;
%     dPCA_results{ds}.per_task_var_per_marg = per_task_var_per_marg;
    
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SUMMARY PLOTS


% 1. MEAN + SD PER-TASK VARIANCE IN THE DPCA MANIFOLD -Pooled across all tasks
% and sessions for each monkey separtely

% get a matrix that has all the values
sptvm = cell2mat(cellfun(@(x) x.per_task_var(:,end), dPCA_results, 'UniformOutput', false )');

figure, hold on
errorbar([1 2],[mean(sptvm(jaco_tasks)) mean(sptvm(jango_tasks))],[std(sptvm(jaco_tasks)) std(sptvm(jango_tasks))],...
    'linestyle','none','linewidth',2,'color','k')
bar([1 2],[mean(sptvm(jaco_tasks)) mean(sptvm(jango_tasks))],'Facecolor',[.6 .6 .6],'linewidth',2)
xlim([0 3]), ylim([0 100])
ylabel('Per-task neural variance expl. (%)')
set(gca,'Tickdir','out'),set(gca,'FontSize',14), box off,
set(gca,'XTick',[1 2],'XTickLabel',{'C','J'})

% 2. MEAN + SD PER-TASK VARIANCE PER MARGINALIZATION IN THE DPCA MANIFOLD