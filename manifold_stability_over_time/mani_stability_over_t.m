%
% Assess manifold stability over time
%
%
%
% Notes:
% - The current version operates on smoothed FRs
%
%

function out = mani_stability_over_t( ds, varargin )


% Default parameters
params                  = struct('mani_dim',        10, ...
                                'nbr_chs',          'min', ...
                                'nbr_iter',         100, ...
                                'task',             'ball', ...
                                'dim_red_method',   'pca', ...
                                'plot_yn',          true);


% read input parameters and replace default values where necessary
param_names         = fieldnames(params);
for p = reshape(varargin,2,[])
    if any(strcmp(p{1},param_names))
        params.(p{1})   = p{2};
    else
        error([p{1} ' is not a recognized parameter name']);
    end
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% A. PREPARE THE DATA

% 0. make sure that ds is 1 x n
if size(ds,1)>size(ds,2)
    ds                  = ds'; 
end

% 1. get the number of neural channels that the user wants
if ischar(params.nbr_chs)
    % get minimum number of channels across datasets
    if strcmp(params.nbr_chs,'min')
        params.nbr_chs  = min( cellfun(@(x) numel(x.neural_chs), ds) );
    end
else
    % check that the specified number of channels is within limits
    if params.nbr_chs > min( cellfun(@(x) numel(x.neural_chs), ds) )
        error('number of channels cannot be greater than the available channels');
    end
end


% 2. find the position of the task in each dataset
aux_task                = cellfun( @(x) strncmp(x.labels,params.task,...
                            length(params.task)), ds, 'UniformOutput', false );
task_pos                = cellfun( @(x) find(x), aux_task )';


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% B. DO PRINCIPAL ANGLES

session_pairs           = nchoosek(1:numel(ds),2);
nbr_session_pairs       = size(session_pairs,1);
pangles                 = cell(1,nbr_session_pairs);
pangles                 = cellfun(@(x) zeros(params.nbr_iter,params.mani_dim), ...
                            pangles, 'UniformOutput', false );


for s = 1:nbr_session_pairs
    for i = 1:params.nbr_iter
        
        % -----------------------------------------------------------------
        % A) Get the data for this iteration
        % get session pairs
        sess1           = session_pairs(s,1);
        sess2           = session_pairs(s,2);
        
        % find position of this task for these sessions
        pos1            = find(task_pos(sess1,:));
        pos2            = find(task_pos(sess2,:));
        
        % get smoothed firing rates --note that column is time
        sfr1            = ds{sess1}.smoothed_FR{pos1}(:,2:end);
        sfr2            = ds{sess2}.smoothed_FR{pos2}(:,2:end);
        
        
        % -----------------------------------------------------------------
        % B) choose a subset of all the channels and do dimensionality
        % reduction
        % B.1. do for dataset 1
        if size(sfr1,2) > params.nbr_chs
            chs1        = datasample( 1:size(sfr1,2), params.nbr_chs, ...
                            'Replace', false );
            chs_disc    = setdiff( 1:size(sfr1,2), chs1 );
            dr1         = dim_reduction( [zeros(size(sfr1,1),1), sfr1], ...
                            params.dim_red_method, chs_disc );
        else
            dr1         = ds{sess1}.dim_red_FR{pos1};
        end
        % B.2. do for dataset 2
        if size(sfr2,2) > params.nbr_chs
            chs2        = datasample( 1:size(sfr2,2), params.nbr_chs, ...
                            'Replace', false );
            chs_disc    = setdiff( 1:size(sfr2,2), chs2 );
            dr2         = dim_reduction( [zeros(size(sfr2,1),1), sfr2], ...
                            params.dim_red_method, chs_disc );
        else
            dr2         = ds{sess2}.dim_red_FR{pos2};
        end
        
        % -----------------------------------------------------------------
        % C) Compute principal angles
        
        pangles{s}(i,:) = principal_angles( dr1.w(:,1:params.mani_dim), ...
                            dr2.w(:,1:params.mani_dim) );
    end
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% C. OUTPUT VARIABLES

out.pangles.raw         = cellfun( @(x) deg2rad(x), pangles, ...
                            'UniformOutput', false );
out.pangles.m           = cellfun( @(x) rad2deg(mean(x,1)), pangles, ...
                            'UniformOutput', false );
out.pangles.sd          = cellfun( @(x) rad2deg(std(x,0,1)), pangles, ...
                            'UniformOutput', false );



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% C. PLOT

if params.plot_yn
    colors              = parula(nbr_session_pairs);

    figure,
    hold on
    for s = 1:nbr_session_pairs
        plot(out.pangles.m{s},'linewidth',2,'color',colors(s,:))
        plot(out.pangles.m{s}+out.pangles.sd{s},'-.','linewidth',1.5,'color',colors(s,:))
        plot(out.pangles.m{s}-out.pangles.sd{s},'-.','linewidth',1.5,'color',colors(s,:))
    end
    xlabel('Manifold dimension'), xlim([0 params.mani_dim+1])
    ylabel('Principal angle (deg)'), ylim([ 0 90])
    set(gca,'TickDir','out','FontSize',14), box off
end
