%
% Plot targets of a current behavior
%
%   function [nbr_targets, target_coord] = plot_targets( binned_data )
%
% Inputs:
%   binned_data         : binned_data struct. It has to contain field
%                           'target' with 'corners' (time, ULx ULy LRx LRy) 
%   task                : task (implemented so far: 'iso', 'spr','wm','ball','mg') 
%   (plot_yn)           : [false] plot the targets
%
% Ouputs:
%   nbr_targets         : number of targets
%   target_coord        : target coordinates (ULx ULy Width Height) for
%                           each target
%
%

function [nbr_targets, target_coord] = get_targets( binned_data, task, varargin )

if nargin == 3
    plot_yn             = varargin{1};
else 
    plot_yn             = false;
end

% find targets
switch task{1}
    case {'wf','iso','spr','iso8','wm'}
        targets         = unique(binned_data.trialtable(:,2:5),'rows');
    case 'ball'
        targets         = [-1 1 1 -1]; % arbitrarily create a target centered in [0 0]
    case {'mg','mg-pt'}
%         % old and new datasets have different trial tables, but we can
%         % distinguish them based on their number of columns. The old trial
%         % table only has 7 columns, and the new one has 12. Unfortunately
%         % in an intermediate time, the target coordinates in the 12-column
%         % trial table were [-1 -1 -1 -1]
%         nbr_cols_mg     = size(binned_data.trialtable,2);
%         switch nbr_cols_mg
%             case 12
        targets = unique(binned_data.trialtable(:,7:10),'rows');
% datasets from Theo don't have target coordinates, it's a
% vector with -1's.
        if targets(1,:) == [-1 -1 -1 -1]
            nbr_targets = numel(find(unique(binned_data.trialtable(:,6))>=0));
            targets = zeros(nbr_targets,4);
            target_h    = 5;
            target_w    = 5;
            for i = 1:nbr_targets
                targets(i,:) = [-1*target_w/2, i*target_h, 1*target_w/2, (i-1)*target_h];
            end
        end
%             case 7
%                 nbr_targets = length(find(unique(binned_data.trialtable(:,4))>=0));
% %                nbr_targets = length(find(unique(binned_data.trialtable(:,6))>=0));
%                 warning('old MG trial table, drawing targets anywhere')
%                 targets = zeros(nbr_targets,4);
%                 for i = 1:nbr_targets
%                     targets(i,:) = [-1 1+i*2 1 -1+i*2];
%                 end
        
end


if ~exist('nbr_targets','var')
    nbr_targets         = size(targets,1);
end


% find bottom left corner (X and Y), width and height for rectangle command
rect_coord              = zeros(nbr_targets,4);
rect_coord(:,1)         = targets(:,1);
rect_coord(:,2)         = min(targets(:,2),targets(:,4));
rect_coord(:,3)         = abs(targets(:,1)-targets(:,3));
rect_coord(:,4)         = abs(targets(:,2)-targets(:,4));

% get rid of targets with width or height zero, if there's any
rect_coord(rect_coord(:,3) == 0,:) = [];
rect_coord(rect_coord(:,4) == 0,:) = [];

% return variables
target_coord            = rect_coord;


% -------------------------------------------------------------------------
% plot
if plot_yn

    colors           	= parula(nbr_targets);
    max_coord         	= max(max(abs(rect_coord)));
    max_y_coord         = max(rect_coord(:,2));
    min_y_coord         = min([min(rect_coord(:,2)), -2]);
    max_height          = max(rect_coord(:,4));
    
    figure
    hold on
    for tg = 1:nbr_targets
        rectangle('Position',rect_coord(tg,:),'Edgecolor',colors(tg,:),...
            'Facecolor',colors(tg,:))    
    end
    rectangle('Position',[-1,-1,2,2],'Edgecolor','k');
    xlim([-max_coord-3, max_coord+3])
    ylim([min_y_coord-max_height, max_y_coord+max_height])
    set(gca,'TickDir','out'),set(gca,'FontSize',14)
    title(['target positions ' task],'FontSize',14);
end