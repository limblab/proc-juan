%
%   comp_nbr                : 1-D vector with the components that will be
%                               compared, or,
%                             2-D vector that established which transformed
%                             dimensions in 'within' will be compared to
%                             which dimensions in 'across'
%   

function [pc_proj_across_tasks, pc_proj_across_tasks_pair2] = transform_and_compare_dim_red_data_all_tasks( dim_red_FR, ...
                    smoothed_FR, labels, neural_chs, comp_nbr, mode, varargin )
    

% -------------------------------------------------------------------------                
% read inputs
if nargin == 7
    if ~iscell(varargin{1})
        plot_yn         = varargin{1};
    else
        % this has to be results.data{task},angles.pair_min_angle,
        % generated with batch_angle_btw_manifolds
        pair_tasks_min_angle = varargin{1};
    end
end
clear varargin;


if ~exist('plot_yn','var')
    plot_yn             = false;
end


% -------------------------------------------------------------------------
% some definitions

% number of tasks
nbr_bdfs                = length(dim_red_FR);
% number of total dimenions in the space (neural chs)
nbr_dims_space          = length(dim_red_FR{1}.eigen);

% possibnle combinations of tasks
comb_bdfs               = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs           = size(comb_bdfs,1);



% -------------------------------------------------------------------------
% compare within manifold and across manifold projections

for i = 1:nbr_comb_bdfs
    
    task_1            	= comb_bdfs(i,1);
    task_2             	= comb_bdfs(i,2);
    
    
    % ---------------------------------------------------------------------
    % choose how the projections will be compared:
    %   1) 'ranking' will compare projections onto eigenvector 'n' in task
    %   1 to projections onto eigenvector 'n' in task 2 (for n = 1:N)
    %   2) 'min_angle' will resort the eigenvectors in one task based on
    %   finding the eigenvector that has the smallest angle with each
    %   eigenvector of the other task
    switch mode
        
        % Define pairs of eigenvectors simply based on its eigenvalue
        % ranking 
        case 'ranking'
            pc_proj_across_tasks(i) = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, task_1, ...
                            task_2, comp_nbr );
        
        % Define pairs of eigenvectors minimizing their angle
        case 'min_angle'
           
            % Reorder the eigenvectors of one of the task based on finding
            % the pairs that define the smallest angle. Note that this
            % operation is not symmetrical, i.e. using task 1 as reference
            % will not yield the same eigenvector ranking as using task 2
            % as reference. In 'comp_neural_manifolds.m' this is solved by
            % doing both and choosing the reference that provides the
            % smallest angle between manifolds.  
            % If the user has passed this information it can be used to
            % resort the eigenvectors to provide maximum alignment.
            % Otherwise it will be arbitrarily chosen
            if exist('pair_tasks_min_angle','var')
                % see what task in the pair we would use as reference for
                % sortin (this would be in column 1
                labels_this_pair = [labels(comb_bdfs(i,1)) labels(comb_bdfs(i,2))];
                % and get the pair that yielded the minimum angle to see if
                % we have to reverse the order
                pair_min_angle = pair_tasks_min_angle{i};
               
                if sum(strcmpi(labels_this_pair,pair_min_angle)) == 2
                    % if the labels match, do not reverse the order
                    disp('ref: task 1')
                    [~, dim_min_angle]  = find_closest_neural_dim( ...
                            dim_red_FR, 1:comp_nbr, '', false );
                        
                    % -----------------------------------------------------
                    % ToDo: delete from here this is to do the other pair
                    % and compare
                    [~, dim_min_angle_pair2]  = find_closest_neural_dim( ...
                            dim_red_FR, 1:comp_nbr, '', true );
                    % ToDo: delete up to here
                    % -----------------------------------------------------
                else
                    % otherwise reverse them
                    disp('ref: task 2')
                    [~, dim_min_angle]  = find_closest_neural_dim( ...
                            dim_red_FR, 1:comp_nbr, '', true );
                        
                    % -----------------------------------------------------
                    % ToDo: delete from here this is to do the other pair
                    % and compare
                    [~, dim_min_angle_pair2]  = find_closest_neural_dim( ...
                            dim_red_FR, 1:comp_nbr, '', false );
                    % ToDo: delete up to here
                    % -----------------------------------------------------
                end
            % this is the case when we don't know if we should use task 1
            % or task 2 as reference
            else
                [~, dim_min_angle]  = find_closest_neural_dim( dim_red_FR, ...
                                1:comp_nbr, '', false );
            end
            
            % format dim_min_angle so it's a 2-D matrix with dimensions in
            % original space, and matching  dimensions in the second space
            [r_eigen, c_eigen] = find(cellfun(@(x) ~isempty(x), dim_min_angle));
            eigenv_pairs    = dim_min_angle{r_eigen,c_eigen};


            % compare within and across projections
            pc_proj_across_tasks(i) = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, task_1, ...
                            task_2, eigenv_pairs, plot_yn );
                        

            % ------------------------------------------------------------ 
            % ToDo: delete from here this is to do the other pair
            % and compare
            
            % this will do the same anlaysis reordering the eigenvectors
            % the opposite way
            
            [r_eigen_2, c_eigen_2] = find(cellfun(@(x) ~isempty(x), dim_min_angle_pair2));
            eigenv_pairs_2  = dim_min_angle_pair2{r_eigen_2,c_eigen_2};
            
            pc_proj_across_tasks_pair2(i) = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, task_1, ...
                            task_2, eigenv_pairs_2, plot_yn );
                        
            % ToDo: delete up to here
            % ------------------------------------------------------------
               
            clear eigenv_pairs
    end
end

