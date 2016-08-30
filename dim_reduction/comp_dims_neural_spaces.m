% 
% Compute the angle between the eigenvectors that define the manifold of
% all pairs of tasks in a BDF struct.
%
%
% Outputs:
%   angle           : cell (nbr_bdfs by nbr_bdfs) with each field being the
%                       angle between eigenvectors. I made this matrix
%                       symmetric although it doesn't need to be.
%   dim_min_angle   : cell (nbr_bdfs by nbr_bdfs) with each field being the
%                       dimension in one task that is closest to each
%                       dimension in the other task. I made this matrix
%                       symmetric although it doesn't need to be.
%   ref_task        : cell (nbr_bdfs by nbr_bdfs) with the reference task
%                       in each pairwise comparison
%
%

function [angle, dim_min_angle, ref_task_label, data] = comp_dims_neural_spaces( ...
    dim_red_FR, dims, labels, varargin )



% PARAMS
plot_all_yn         = false;


% -------------------------------------------------------------------------
% Read inputs

if nargin == 4
    angle_orth      = varargin{1};
end


% if the angle at which angles become non-orthogonal has not been passed,
% compute it according to the eq. in the Suppl. Info of (Kobak et al.,
% eLife 2016)
if ~exist('angle_orth','var')
    angle_orth      = rad2deg(acos(3.3/sqrt( length(dim_red_FR{1}.eigen) )));
end
if isempty(angle_orth)
    angle_orth      = rad2deg(acos(3.3/sqrt( length(dim_red_FR{1}.eigen) )));
end


% -------------------------------------------------------------------------
% Definitions

nbr_bdfs            = length(dim_red_FR);

% matrix with all possible pairs of tasks
comb_bdfs           = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs       = size(comb_bdfs,1);


% define cells for storing the results
data.angle          = cell(nbr_bdfs);
data.dim_min_angle  = cell(nbr_bdfs);
angle               = cell(nbr_bdfs);
dim_min_angle       = cell(nbr_bdfs);
ref_task_label      = cell(nbr_bdfs);


% -------------------------------------------------------------------------
% Compute angles between all pairs of manifolds, using both tasks in each
% pair as reference to align them


% do for all pairs of tasks
for p = 1:nbr_comb_bdfs
    
    % Order eigenvectors and compute angle using task 1 in the pair of
    % tasks as reference 
    [ang, dim_m_ang] = find_closest_neural_dim( dim_red_FR, dims );

    % Do the same using task 2 as reference
    [ang_rev, dim_m_ang_rev] = find_closest_neural_dim( dim_red_FR, dims, [], true );
    
    
    % store angles in return cells 
    data.angle{comb_bdfs(p,1),comb_bdfs(p,2)} = ang{comb_bdfs(p,1),comb_bdfs(p,2)};
    data.angle{comb_bdfs(p,2),comb_bdfs(p,1)} = ang_rev{comb_bdfs(p,2),comb_bdfs(p,1)};
    
    % store dimensions that give the smallest angle
    data.dim_min_angle{comb_bdfs(p,1),comb_bdfs(p,2)} = dim_m_ang{comb_bdfs(p,1),comb_bdfs(p,2)};
    data.dim_min_angle{comb_bdfs(p,2),comb_bdfs(p,1)} = dim_m_ang_rev{comb_bdfs(p,2),comb_bdfs(p,1)};
end



% -------------------------------------------------------------------------
% Decide which is the best task to use as reference, by looking at how many
% pairs of eigenvectors are not orthogonal -- the one that gives the
% smallest number will be used as reference

for p = 1:nbr_comb_bdfs
    
    nbr_orth_eigenv_1   = numel(find(rad2deg(data.angle{comb_bdfs(p,1),comb_bdfs(p,2)}) > angle_orth ));
    nbr_orth_eigenv_2   = numel(find(rad2deg(data.angle{comb_bdfs(p,2),comb_bdfs(p,1)}) > angle_orth ));
    
    
    % -----------------------------
    % choose the task with the smallest number of eigenvectors above the
    % "orthogonality" threshold. If they are the same, choose the one that
    % goes above the threshold last
    if nbr_orth_eigenv_1 < nbr_orth_eigenv_2
        ref_task        = 2;
        non_ref_task    = 1;
    elseif nbr_orth_eigenv_1 > nbr_orth_eigenv_2
        ref_task        = 1;
        non_ref_task    = 2;
    else % same number of non-orthogonal eigenvectors
        if find(rad2deg(data.angle{comb_bdfs(p,1),comb_bdfs(p,2)}) > angle_orth, 1) > ...
                find(rad2deg(data.angle{comb_bdfs(p,2),comb_bdfs(p,1)}) > angle_orth, 1 )
            ref_task    = 1;
            non_ref_task = 2;
        else
            ref_task    = 2;
            non_ref_task =1;
        end
    end
    
    
    % -----------------------------
    % Store results
    
    angle{comb_bdfs(p,1),comb_bdfs(p,2)} = data.angle{comb_bdfs(p,ref_task),comb_bdfs(p,non_ref_task)};
    dim_min_angle{comb_bdfs(p,1),comb_bdfs(p,2)} = data.dim_min_angle{comb_bdfs(p,ref_task),comb_bdfs(p,non_ref_task)};
    ref_task_label{comb_bdfs(p,1),comb_bdfs(p,2)} = labels{comb_bdfs(p,ref_task)};
        
    % make matrices diagonal
    angle{comb_bdfs(p,2),comb_bdfs(p,1)} = angle{comb_bdfs(p,1),comb_bdfs(p,2)};
    dim_min_angle{comb_bdfs(p,2),comb_bdfs(p,1)} = dim_min_angle{comb_bdfs(p,1),comb_bdfs(p,2)};
    ref_task_label{comb_bdfs(p,2),comb_bdfs(p,1)} = ref_task_label{comb_bdfs(p,1),comb_bdfs(p,2)};
end



% -------------------------------------------------------------------------
% Plot

% plots for all tasks
if plot_all_yn
    for p = 1:nbr_comb_bdfs
        figure,hold on
        plot(rad2deg(data.angle{comb_bdfs(p,1),comb_bdfs(p,2)}),'linewidth',2,'marker','d')
        plot(rad2deg(data.angle{comb_bdfs(p,2),comb_bdfs(p,1)}),'linewidth',2,'marker','d','color','r')
        plot([0 dims(end)+1],[angle_orth angle_orth],':','color',[.6 .6 .6],'linewidth',2)
        plot(rad2deg(angle{comb_bdfs(p,1),comb_bdfs(p,2)}),':k','linewidth',2,'marker','d')
        legend(['ref. ' labels{comb_bdfs(p,1)}],['ref. ' labels{comb_bdfs(p,2)}],'orth (P>0.001)','chosen','Location','SouthEast'),legend boxoff
        xlabel('dimension'), ylabel('angle')
        set(gca,'TickDir','out'), set(gca,'FontSize',16)
        xlim([0 dims(end)+1]), ylim([0 90])

        figure,hold on
        plot(data.dim_min_angle{comb_bdfs(p,1),comb_bdfs(p,2)}(:,2),...
            data.dim_min_angle{comb_bdfs(p,2),comb_bdfs(p,1)}(:,2),'k','linewidth',2,'marker','d')
        plot([0 50],[0 50],'color',[.5 .5 .5])
        xlabel(['ref. ' labels{comb_bdfs(p,1)}]),ylabel(['ref. ' labels{comb_bdfs(p,2)}])
        legend('closest eigenv.','1:1','Location','SouthEast')
        set(gca,'TickDir','out'), set(gca,'FontSize',16)
    end
end

% summary plot
figure,hold on
aux_lgnd            = cell(1,nbr_comb_bdfs+1);
aux_colors          = distinguishable_colors(nbr_comb_bdfs);
for p = 1:nbr_comb_bdfs
    plot(rad2deg(angle{comb_bdfs(p,1),comb_bdfs(p,2)}),'color',aux_colors(p,:),'linewidth',2,'marker','d')
    aux_lgnd{p}     = [ref_task_label{comb_bdfs(p,1),comb_bdfs(p,2)} ' vs. ' ...
        labels{comb_bdfs(p,find(strcmpi(labels(comb_bdfs(p,:)),ref_task_label(comb_bdfs(p,1),comb_bdfs(p,2)))==0))}];
end
plot([0 dims(end)+1],[angle_orth angle_orth],':','color',[.6 .6 .6],'linewidth',2)
aux_lgnd{end}       = 'orth (P>0.001)';
legend(aux_lgnd,'Location','SouthEast'),legend boxoff
xlabel('dimension'), ylabel('angle')
set(gca,'TickDir','out'), set(gca,'FontSize',16)
xlim([0 dims(end)+1]), ylim([0 90])

end