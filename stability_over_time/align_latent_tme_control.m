%
% Function that uses TME to generate surrogate datasets and see how well we
% can align surrogate neurons that preserve the tuning across target
% directions and over time
%

function results = align_latent_tme_control( td, params )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES

signals             = [];
xval_yn             = false;
n_folds             = 5;
method              = 'cca'; % 'cca' or 'procrustes'
mani_dims           = 1:10;
surrogate_type      = 'surrogate-C';
n_surrogates        = 1000;
mani_dims           = [];

if nargin > 1, assignParams(who,params); end % overwrite defaults

if ~strcmp(method,'cca'); error('Procrustes not implemented yet! only PCA'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all pairwise comparisons of sessions

sessions            = unique({td.date});
n_sessions          = length(sessions);

comb_sessions       = nchoosek(1:n_sessions,2);
n_comb_sessions     = size(comb_sessions,1);


% overwrite signals for the spikes
spiking             = [signals(1:end-4) '_spikes'];

rng('shuffle', 'twister') % randomize the seed for TME


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align using CCA or Procrustes


results.cc_mn       = zeros(n_comb_sessions,numel(mani_dims));
results.cc_sd       = zeros(n_comb_sessions,numel(mani_dims));
results.prc99       = zeros(n_comb_sessions,numel(mani_dims));

for c = 1:n_comb_sessions

    
    s1              = comb_sessions(c,1);
    s2              = comb_sessions(c,2);
    
    [~,td1]         = getTDidx( td, {'date', sessions{s1}} );
    [~,td2]         = getTDidx( td, {'date', sessions{s2}} );
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Prepare the data
    
    % Trial average the data
    td1avg          = trialAverage( td1, 'target_direction' );
    td2avg          = trialAverage( td2, 'target_direction' );
    
    n_bins          = size(td1avg(1).pos,1);
    
    % Assemble tensor:  time x unit x condition
    fr1             = zeros( n_bins, size(td1avg(1).(spiking),2), length(td1avg) );
    fr2             = zeros( n_bins, size(td2avg(1).(spiking),2), length(td2avg) );
    
    for t = 1:length(td1avg)
       fr1(:,:,t)   = td1avg(t).(spiking);
       fr2(:,:,t)   = td2avg(t).(spiking);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Prepare for TME
    
    % Quantify primary features of the original data
    [targetSigmaT1, targetSigmaN1, targetSigmaC1, M1] = extractFeatures_J(fr1);
    [targetSigmaT2, targetSigmaN2, targetSigmaC2, M2] = extractFeatures_J(fr2);
    
    
    % define surrogate params
    switch surrogate_type
        case 'surrogate-TC'
            params1.margCov{1} = targetSigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = targetSigmaC1;
            params1.meanTensor = M1.TC;

            params2.margCov{1} = targetSigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = targetSigmaC2;
            params2.meanTensor = M2.TC;
        case 'surrogate-C'
            params1.margCov{1} = [];
            params1.margCov{2} = [];
            params1.margCov{3} = targetSigmaC1;
            params1.meanTensor = M1.C;

            params2.margCov{1} = [];
            params2.margCov{2} = [];
            params2.margCov{3} = targetSigmaC2;
            params2.meanTensor = M2.C;
        case 'surrogate-T'
            params1.margCov{1} = targetSigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = [];
            params1.meanTensor = M1.T;

            params2.margCov{1} = targetSigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = [];
            params2.meanTensor = M2.T;
        otherwise
            error('Need to code these TME parameters')
    end

    
    % Fit the maximum entropy distribution
    maxEntropy1     = fitMaxEntropy(params1);
    maxEntropy2     = fitMaxEntropy(params2);

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Run TME

    disp(['Generating surrogate distribution for session comparison ' num2str(c) ' of ' num2str(n_comb_sessions)]);
    
    surr_align      = zeros(n_surrogates,numel(mani_dims));

    for s = 1:n_surrogates
        
        % generate surrogate datasets
        surr_tensor1 = sampleTME(maxEntropy1);
        surr_tensor2 = sampleTME(maxEntropy2);

        % Do PCA of the surrogate datasets
        surr_tensor1 = permute(surr_tensor1,[2 1 3]);
        surr_tensor2 = permute(surr_tensor2,[2 1 3]);
        
        surr_fr1    = reshape(surr_tensor1, size(surr_tensor1,1), [])';
        surr_fr2    = reshape(surr_tensor2, size(surr_tensor2,1), [])';
        
        [~,sc1]     = pca(surr_fr1);
        [~,sc2]     = pca(surr_fr2);

        [~,~,surr_align(s,:)] = canoncorr( sc1(:,mani_dims), sc2(:,mani_dims) );
    end
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DO THE SAME FOR TRIAL AVERAGED DATA

    sp1             = cat(1,td1avg.(spiking));
    sp2             = cat(1,td2avg.(spiking));

    [~,sc1]         = pca(sp1);
    [~,sc2]         = pca(sp2);

    [~,~,actual_cc(c,:)] = canoncorr(sc1(:,mani_dims),sc2(:,mani_dims));
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store results
    
    results.cc_mn(c,:)  = mean(surr_align,1);
    results.cc_sd(c,:)  = std(surr_align,0,1);
    results.prc99(c,:)  = prctile(surr_align,99);
    results.data{c}     = surr_align;
end


% also return actual CCs
results.cc_actual_trialavg = actual_cc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of days between sessions, to summarize the results

diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end

results.diff_days       = diff_days;
results.comb_sessions   = comb_sessions; 





ccs_plot = 1:4;

mn_actual = mean(actual_cc(:,ccs_plot),2);
sd_actual = std(actual_cc(:,ccs_plot),0,2);

mn_surr = mean(results.cc_mn(:,ccs_plot),2);
sd_surr = std(results.cc_mn(:,ccs_plot),0,2);


figure,hold on
errorbar(diff_days,mn_actual,sd_actual,'.k','markersize',32,'linestyle','none')
errorbar(diff_days,mn_surr,sd_surr,'.','color',[1 .6 0],'markersize',32,'linestyle','none')
set(gcf, 'color', [1 1 1])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days between sessions'), 
ylabel(['Similarity latent activity (top ' num2str(numel(ccs_plot)) ' modes) (r)'])
xlim([0 max(diff_days)+1]),ylim([0 1.1])
title(surrogate_type)
legend('actual','surr','Location','SouthWest');legend boxoff
