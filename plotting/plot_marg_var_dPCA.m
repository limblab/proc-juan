%
% Bar plot summarizing the marginalized variance per condition for one or
% several dPCA datasets (in a cell array)
%

function plot_marg_var_dPCA( dPCA_results, varargin ) 

if nargin == 1
    marg_labels         = {'task','target','time','task/target'};
else
    marg_labels         = varargin{1};
end
   

nbr_sessions            = length(dPCA_results);


% if there are multiple sessions
if nbr_sessions > 1
    
    % get the nbr of marginalizations and check that we have the right nbr
    % of labels
    nbr_margs           = length(dPCA_results{1}.expl_var.totalMarginalizedVar);
    if length(marg_labels) ~= nbr_margs
        error('nbr. labels does not match nbr. of marginalizations!');
    end
    
    % fill matrices with expl var per marginalization (in %)
    expl_var            = zeros(nbr_sessions,nbr_margs);
    
    for i = 1:nbr_margs
        expl_var(:,i)   = cell2mat( cellfun(@(x) x.expl_var.totalMarginalizedVar(i)/x.expl_var.totalVar, ...
                            dPCA_results, 'UniformOutput', false)' )*100;
    end
    

% if there's only one session
else

    % get the nbr of marginalizations and check that we have the right nbr
    % of labels
    nbr_margs           = length(dPCA_results.expl_var.totalMarginalizedVar);
    
    if length(marg_labels) ~= nbr_margs
        error('nbr. labels does not match nbr. of marginalizations!');
    end
    
    % fill matrices with expl var per marginalization (in %)
    expl_var            = zeros(1,nbr_margs);
    
    for i = 1:nbr_margs
        expl_var(i)     = dPCA_results.expl_var.totalMarginalizedVar(i)/dPCA_results.expl_var.totalVar*100;
    end
end


% PLOT
figure, hold on
errorbar(1:nbr_margs,mean(expl_var,1),std(expl_var,0,1),'k','linewidth',2,'linestyle','none')
bar(1:nbr_margs,mean(expl_var,1))
set(gca,'TickDir','out'), set(gca,'FontSize',16)
set(gca,'XTick',1:nbr_margs)
set(gca,'XTickLabel',marg_labels)
ylabel('explained variance (%)')