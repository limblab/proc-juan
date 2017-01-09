%
% Plot targets of a current behavior
%
%   function [nbr_targets, target_coord] = get_targets( binned_data, task )
%
% Inputs:
%   binned_data         : binned_data struct. It has to contain field
%                           'target' with 'corners' (time, ULx ULy LRx LRy) 
%   task                : task (implemented so far: 'iso', 'spr','wm','ball','mg') 
%   (plot_yn)           : [false] plot the targets
%
% Ouputs:
%   nbr_targets         : number of targets
%   target_coord        : 2D matrix with each row being the coordinates of
%                           each of the targets (ULx ULy LRx LRy)
%   target_nbr          : 1D matrix with the sequence in which the targets
%                           appeared. Targets are named 1,2, ..., n, where
%                           n is the number of rows in target_coord
%   target_coord_mtlb   : 2D matrix with target coordinates (ULx ULy Width
%                           Height) for each target --this is the format
%                           matlab likes for plotting them with rectangle
%
%

function [nbr_targets, target_coord, target_nbr, varargout] = ...
            get_targets( binned_data, task, varargin )

if nargin == 3
    plot_yn             = varargin{1};
else 
    plot_yn             = false;
end

% for compatibility with different wrapper functions
if iscell(task)
    task_name           = task{1};
else
    task_name           = task;
end


% find targets
switch task_name
    case {'wf','iso','spr','iso8','wm'}
        targets         = unique(binned_data.trialtable(:,2:5),'rows');
    case 'ball'
        targets         = [-1 1 1 -1]; % arbitrarily create a target centered in [0 0]
    case {'mg','mg-pt'}
        % old and new datasets have different trial tables, but we can
        % distinguish them based on their number of columns. The old trial
        % table only has 7 columns, and the new one has 12. Unfortunately
        % in an intermediate time, the target coordinates in the 12-column
        % trial table were always [-1 -1 -1 -1]
%         nbr_cols_mg     = size(binned_data.trialtable,2);
%         switch nbr_cols_mg
%             case 12
        targets = unique(binned_data.trialtable(:,7:10),'rows');
        % datasets from Theo don't have target coordinates, it's a vector
        % with -1's. The code below just converts them to some coordinates
        % that make sense for plotting
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


% -------------------------------------------------------------------------
% get target coords in Matlab's plotting format

% find bottom left corner (X and Y), width and height for rectangle command
rect_coord              = zeros(nbr_targets,4);
rect_coord(:,1)         = targets(:,1);
rect_coord(:,2)         = min(targets(:,2),targets(:,4));
rect_coord(:,3)         = abs(targets(:,1)-targets(:,3));
rect_coord(:,4)         = abs(targets(:,2)-targets(:,4));

% get rid of targets with width or height zero, if there's any
rect_coord(rect_coord(:,3) == 0,:) = [];
rect_coord(rect_coord(:,4) == 0,:) = [];



% -------------------------------------------------------------------------
% create vector with target ID for each trial

nbr_trials              = size(binned_data.trialtable,1);

% FIX for some old wrist-flexion task: sometimes, the go time == -1 instead
% of the real time, in which case we won't count that trial, and we will
% remove it from the trial table
if find( strcmpi(task,{'wf','iso','spr','iso8','wm'}) )
    wrong_time_trials   = find(binned_data.trialtable(:,6)==-1);
    nbr_trials          = nbr_trials - numel(wrong_time_trials);
    binned_data.trialtable(wrong_time_trials,:) = [];
end


target_nbr              = zeros(1,nbr_trials);

switch task_name
    case {'wf','iso','spr','iso8','wm'}
        
        target_indxs    = arrayfun(@(x,y,z,w) find(binned_data.trialtable(:,2)==x & binned_data.trialtable(:,3)==y & binned_data.trialtable(:,4)==z & binned_data.trialtable(:,5)==w), ...
            targets(:,1), targets(:,2), targets(:,3), targets(:,4),'UniformOutput',false);
        % ToDo: double check this DIRTY FIX
        % --it seems that sometimes there is one target that is mistakenly
        % labelled, and it's only one. The code will get rid of it
        indx_discard    = cellfun( @(x) find(length(x)==1), target_indxs, 'UniformOutput', false ); 
        target_indxs(find(cellfun(@(x) ~isempty(x), indx_discard))) = [];
        
        % update the number of targets
        nbr_targets     = length(target_indxs);
        
        % fill a vector of target numbers with a scalar that identifies the
        % target (>= 0)
        target_nbr      = nan(1,sum(cellfun(@(x) length(x), target_indxs)));
        for i = 1:length(target_indxs)
            target_nbr(target_indxs{i}) = i-1;
        end
        
    case 'ball'
        target_nbr      = zeros(1,nbr_trials);
    case {'mg','mg-pt'}
        target_nbr      = binned_data.trialtable(:,6);
end


% -------------------------------------------------------------------------
% return variables
if exist('indx_discard','var')
    target_coord        = targets(setdiff(1:size(targets,1),...
                            find(cellfun(@(x) ~isempty(x), indx_discard))),:);
else
    target_coord        = targets;
end
if nargout == 4
    varargout{1}        = rect_coord;
end


% -------------------------------------------------------------------------
% plot
if plot_yn

    %colors           	= distinguishable_colors(nbr_targets);
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