%
% Shuffle elements in covariance matrix 
%

% Params
n_dims = 12;

% dataset and task
ds = 11;
tsk = 1;

% number of repetitions
n_sh = 10000;

% compute covariance matrix of the de-meaned neural data
X = datasets{ds}.stdata{tsk}.target{end}.neural_data.conc_smoothed_fr;
X = X - repmat(mean(X,1),size(X,1),1);

C = cov(X);


% -------------------------------------------------------------------------
% shuffle the off-diagonal elements in the covariance matrix 


% preallocate matrix for randomly generated angles
theta = zeros(n_sh,n_dims);

% % get idx off-diagonal elements
% i_od = find(ones(size(C,1))-eye(size(C,1)));

% get idx of upper diagonal elements
i_ud = find(triu(C,1));

for s = 1:n_sh
    
%     % indexes to shuffle off-diagonal elements
%     i_sh = i_od(randperm(length(i_od)));
% 
%     % get resultant shuffled covariance matrix
%     C_sh = C;
%     for j = 1:length(i_od)
%         C_sh(i_od(j)) = C(i_sh(j));
%     end

    % indexes to shuffle the upper diagonal elements
    i_sh = i_ud(randperm(length(i_ud)));
    
    % get the shuffled covariance matrix
    C_sh = eye(size(C,1)).*C; % keep same diagonal elements
    for j = 1:length(i_ud)
        C_sh(i_ud(j)) = C(i_sh(j));
    end
    C_sh = C_sh + triu(C_sh,1)';
    
    
    % Do SVD of the original and shuffled covariances matrices to "do PCA"
    [U, S, V] = svd(C);
    [U_sh, S_sh, V_sh] = svd(C_sh);
    
%     figure,hold on
%     plot(diag(S)/sum(diag(S)),'b','linewidth',2), plot(diag(S_sh)/sum(diag(S_sh)),'r','linewidth',2)
%     set(gca,'TickDir','out','FontSize',14), box off
%     xlabel('Principal component'),ylabel('Variance explained')
%     legend('Original','Shuffled Cov'),legend boxoff

    % compute principal angles for this pair
    theta(s,:) = principal_angles(V(:,1:n_dims),V_sh(:,1:n_dims));

end



figure,
plot(rad2deg(theta'),'color',[.6 .6 .6]),
hold on
plot(rad2deg(prctile(theta,1)),'k','linewidth',2)
ylim([0 90]),xlim([0 n_dims])
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Neural mode'),ylabel('Principal angle (deg)')

% plot real principal angles
W1 = datasets{ds}.dim_red_FR{tsk}.w(:,1:n_dims);
W2 = datasets{ds}.dim_red_FR{tsk+1}.w(:,1:n_dims);
plot(rad2deg(principal_angles(W1,W2)),'r','linewidth',2)