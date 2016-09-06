

function proj_results = batch_cross_manifold_proj( varargin )



% -------------------------------------------------------------------------
% read input params

if nargin >=1
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end

if nargin >= 2
    angle_results       = varargin{2};
end

if nargin == 3
    params              = varargin{3};
else
    params              = batch_cross_manifold_proj_params_defaults;
end

clear varargin;


% -------------------------------------------------------------------------
% load data

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end


% run 'angle' analysis as we need it for comparing the projections
if ~exist('angle_results','var')
    angle_results       = batch_angle_btw_neural_manifolds( datasets );
end


% -------------------------------------------------------------------------
% get info about monkeys, tasks, and who did what

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% assess how similat the within-manifold and across-manifold projections
% are

for i = 1:meta_info.nbr_monkeys
     
    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
        
        % what dataset (monkey and session) are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);

        % compare projections onto each of the dimensions of the within
        % task and the across task manifold in two cases:
        %
        %   1) ranking the eigenvectors in the way that yielded the highest
        %   dimensional, non-orthogonal manifolds ('pc_projs')
        %   2) using the alternative eigenvector ranking ('pc_projs_pair2')
        %
        [pc_projs, pc_projs_pair2] = transform_and_compare_dim_red_data_all_tasks( ...
                            datasets{dtst}.dim_red_FR, ...
                            datasets{dtst}.smoothed_FR, ...
                            datasets{dtst}.labels, ...
                            datasets{dtst}.neural_chs, ...
                            params.dim_manifold, ...
                            'min_angle', ...
                            angle_results.data{dtst}.angles_manifolds.pair_min_angle );
                        
        % store all results
        proj_comp{dtst} = pc_projs; 
        proj_comp_p2{dtst} = pc_projs_pair2;
        
        
        % plot r and R2 and corresponding eigenvectors for both
        if params.plot_p_pair
            
            for j = 1:length(pc_projs)
                ttl = datasets{dtst}.monkey; 
                ttl = [ttl ' ' char(angle_results.data{dtst}.angles_manifolds.pair_min_angle{j}(1)) ...
                    ' - ' char(angle_results.data{dtst}.angles_manifolds.pair_min_angle{j}(2))];
                figure('units','normalized','outerposition',[1/8 1/4 3/4 1/2]),
                subplot(131),hold on, 
                plot(pc_projs(j).r,'linewidth',3,'color','b'),plot(pc_projs_pair2(j).r,'r','linewidth',3)
                ylabel('r'),xlabel('principal axis'),
    %             plot(find(pc_projs(j).P_r<0.001),ones(1,length(find(pc_projs(j).P_r<0.001))),'linestyle','none','marker','.','markersize',10,'color','b')
    %             plot(find(pc_projs_pair2(j).P_r<0.001),1.05*ones(1,length(find(pc_projs_pair2(j).P_r<0.001))),'linestyle','none','marker','.','markersize',10,'color','r')
                plot(max(pc_projs(j).r_shuffled),'linewidth',3,'color','c','linestyle',':')
                plot(max(pc_projs_pair2(j).r_shuffled),'r','linewidth',3,'color',[1 .6 0],'linestyle',':')
                plot(max(pc_projs(j).r_other_projs),'linewidth',3,'color','c','linestyle','-.')
                plot(max(pc_projs_pair2(j).r_other_projs),'linewidth',3,'color',[1 .6 0],'linestyle','-.')
                legend('min angle','max angle','shuf min angle','shuf max angle',...
                    'max other PC','max other PC','Location','NorthEast'), legend boxoff
                ylim([0 1])
                set(gca,'Tickdir','out'),set(gca,'FontSize',14)
                title(ttl)

                subplot(132),hold on, plot(pc_projs(j).R2,'linewidth',3),plot(pc_projs_pair2(j).R2,'r','linewidth',3)
                plot(max(pc_projs(j).R2_shuffled),'linewidth',3,'color','c','linestyle',':')
                plot(max(pc_projs_pair2(j).R2_shuffled),'r','linewidth',3,'color',[1 .6 0],'linestyle',':')
                plot(max(pc_projs(j).R2_other_projs),'linewidth',3,'color','c','linestyle','-.')
                plot(max(pc_projs_pair2(j).R2_other_projs),'linewidth',3,'color',[1 .6 0],'linestyle','-.')
                set(gca,'Tickdir','out'),set(gca,'FontSize',14)
                ylim([0 1]),ylabel('R^2')

                subplot(133),hold on
                plot([0 30],[0 30],'color',[.6 .6 .6],'linewidth',3)
                plot(pc_projs(j).comp_nbr(:,2),pc_projs_pair2(j).comp_nbr(:,2),'.k','markersize',20),
                xlabel('matched eigenv min angle'),ylabel('matched eigenv max angle')
                set(gca,'Tickdir','out'),set(gca,'FontSize',14)
            end
        end
    end
end


% -------------------------------------------------------------------------
% Return results

for d = 1:length(datasets)
    proj_results{d}.min_manifold_angle_pair = proj_comp{d};
    proj_results{d}.max_manifold_angle_pair = proj_comp_p2{d};
end





% -------------------------------------------------------------------------
% Some checks to see if the dimension that defines the smallest angle also
% gives the maximum cross-task projection similarity

for d = 1:length(datasets)
    
    for p = 1:length(proj_results{d}.min_manifold_angle_pair)
    
        % retrieve correlation
        corr_min_angle      = proj_results{d}.min_manifold_angle_pair(p).r;
        corr_max_angle      = proj_results{d}.max_manifold_angle_pair(p).r;
        
        % find maximum correlation when projecting on other dimensions
        [max_r_other_min_mani_angle, eigenv_max_r_other_min_mani_angle] = ...
                    max(proj_results{d}.min_manifold_angle_pair(p).r_other_projs);
        [max_r_other_max_mani_angle, eigenv_max_r_other_max_mani_angle] = ...
                    max(proj_results{d}.max_manifold_angle_pair(p).r_other_projs);
                
        % find maximum correlation across all these 4 conditions, for each
        % dimension
        aux_r = [corr_min_angle; corr_max_angle; max_r_other_min_mani_angle; max_r_other_max_mani_angle];
        [max_r, pair_max_r] = max(aux_r,[],1);
        
        % plot correlations
        figure,imagesc(aux_r), colorbar
        set(gca,'Tickdir','out'),set(gca,'FontSize',14)
        set(gca,'YTick',1:4), set(gca,'YTickLabelRotation',45)
%        set(gca,'YTickLabel',{'other max r max','other max r min','max mani ang','min mani ang'})
        set(gca,'YTickLabel',{'min mani ang','max mani ang','other max r min','other max r max'})
        hold on, plot(pair_max_r,'linestyle','none','color','w','marker','.','markersize',24)
        title([datasets{d}.monkey ' - ' angle_results.data{d}.angles_manifolds.pair_min_angle{p}{1} ...
            ' vs ' angle_results.data{d}.angles_manifolds.pair_min_angle{p}{2}]);
        
        
%         % compare this max correlation to the correlation onto the direction that
%         % defines the smallest angle
%         indx_higher_r_than_corr_eigenv_min_mani = find( max_r_other_min_mani_angle > ...
%             proj_results{d}.min_manifold_angle_pair(p).r );
%         indx_higher_r_than_corr_eigenv_max_mani = find( max_r_other_max_mani_angle > ...
%                     proj_results{d}.max_manifold_angle_pair(p).r );

%         figure,
%         subplot(121),hold on
%         plot(proj_results{d}.min_manifold_angle_pair(p).r,'linewidth',2,'color','b')
%         plot(proj_results{d}.max_manifold_angle_pair(p).r,'linewidth',2,'color','k')
%         plot(max_r_other_min_mani_angle,'linewidth',2,'color','c')
%         plot(max_r_other_max_mani_angle,'linewidth',2,'color',[.6 .6 .6])
%         legend('min mani angle','max mani angle','opt other min mani','opt other max mani')
%         legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14)
%         xlim([0 10])
%         title([meta_info.tasks_per_session{d}{proj_results{d}.min_manifold_angle_pair(p).within_task}, ...
%             ' vs. ' meta_info.tasks_per_session{d}{proj_results{d}.min_manifold_angle_pair(p).across_task}]);
%         
%         subplot(122),hold on
%         plot(proj_results{d}.min_manifold_angle_pair(p).comp_nbr(:,2),'linewidth',2,'color','b')
%         plot(proj_results{d}.max_manifold_angle_pair(p).comp_nbr(:,2),'linewidth',2,'color','k')
%         plot(eigenv_max_r_other_min_mani_angle,'linewidth',2,'color','c')
%         plot(eigenv_max_r_other_max_mani_angle,'linewidth',2,'color',[.6 .6 .6])
%         legend('min mani angle','max mani angle','opt other min mani','opt other max mani')
%         legend boxoff
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14)
%         xlim([0 10])
    end
end
    



% % -------------------------------------------------------------------------
% % Plot to compare the dimension that gives the smallest angle between manifolds
% % and the dimensions that gives the most correlated cross-task projection
% for d = 1:length(datasets)
%     for p = 1:length(proj_results{d}.min_manifold_angle_pair)
%         figure,hold on
%         plot(proj_results{d}.min_manifold_angle_pair(p).comp_nbr(:,2),'k')
%         plot(proj_results{d}.max_manifold_angle_pair(p).comp_nbr(:,2),'r')
% 
%         % ToDo..
%         
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14)
%     end
% end