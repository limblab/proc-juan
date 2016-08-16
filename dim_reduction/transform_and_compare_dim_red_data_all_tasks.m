%
%   comp_nbr                : 1-D vector with the components that will be
%                               compared, or,
%                             2-D vector that established which transformed
%                             dimensions in 'within' will be compared to
%                             which dimensions in 'across'
%   

function pc_proj_across_tasks = transform_and_compare_dim_red_data_all_tasks( dim_red_FR, ...
                    smoothed_FR, labels, neural_chs, comp_nbr, mode, varargin )
    
                
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
                
% number of tasks
nbr_bdfs                = length(dim_red_FR);
% number of total dimenions in the space (neural chs)
nbr_dims_space          = length(dim_red_FR{1}.eigen);

% possibnle combinations of tasks
comb_bdfs               = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs           = size(comb_bdfs,1);


for i = 1:nbr_comb_bdfs
    
    bdf_1               = comb_bdfs(i,1);
    bdf_2               = comb_bdfs(i,2);
    
    switch mode
        
        case 'ranking'
            pc_proj_across_tasks(i) = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, bdf_1, ...
                            bdf_2, comp_nbr );
             
        case 'min_angle'
            
            
            % ToDo: legecy code --delete
%             [~, dim_min_angle_old] = find_closest_hyperplane( dim_red_FR{bdf_1}.w, ...
%                             dim_red_FR{bdf_2}.w, 1:nbr_dims_space );
            % ToDo: delete until here
            
            
            % Reorder the eigenvectors of one of the task based on finding
            % the pairs that define the smallest angle. Note that this
            % operation is not symmetrical, i.e. using task 1 as reference
            % will not yield the same eigenvector ranking as using task 2
            % as reference. In comp_neural_spaces_fcn_dim_finding_closest
            % this is solved by doing both and choosing the reference that
            % provides the smallest angle between manifolds. 
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
                    disp('same order')
                    [~, dim_min_angle]  = find_closest_neural_hyperplane_all( ...
                            dim_red_FR, 1:comp_nbr, '', false );
                else
                    % otherwise reverse them
                    disp('reversing order')
                    [~, dim_min_angle]  = find_closest_neural_hyperplane_all( ...
                            dim_red_FR, 1:comp_nbr, '', true );
                end
            % this is the case when we don't know if we should use task 1
            % or task 2 as reference
            else
                [~, dim_min_angle]  = find_closest_neural_hyperplane_all( dim_red_FR, ...
                                1:comp_nbr, '', false );
            end
            
            % format dim_min_angle so it's a 2-D matrix with dimensions in
            % original space, and matching  dimensions in the second space
            [r_eigen, c_eigen] = find(cellfun(@(x) ~isempty(x), dim_min_angle));
            eigenv_pairs    = dim_min_angle{r_eigen,c_eigen};

            
% ToDo: legecy code --delete
%             % eigenv_pairs                        = [1:nbr_dims_space;dim_min_angle_old];
%             eigenv_pairs    = dim_min_angle; 
%             
%             % truncate to the number of dimensions we want
%             eigenv_pairs(:,comp_nbr+1:end)      = [];
% ToDo: delete up to here


            % compare within and across projections
            pc_proj_across_tasks(i) = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, bdf_1, ...
                            bdf_2, eigenv_pairs, plot_yn );
               
            clear eigenv_pairs
    end
end

