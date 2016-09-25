%
% Plot single trial data struct
%
%   function plot_single_trial_data( single_trial_data, var, params )
%
%
% Inputs (opt)          : [default]
%   single_trial_data   : single_trial_data cell struct
%   var                 : variable to plot ('neuron', 'emg', 'neural PC',
%                           'muscle synergy', 'pos', 'vel', 'force') 
%   param               : what neurons/emgs/pcs/synergies/pos signals/vel
%                           signals/force signals to plot
%   (target)            : ['each'] target number: 'each') one plot per
%                           target, showing the data from all the trials,
%                           as well as the mean (+/-SD); scalar 1:n) target
%                           number, same format as 'each'; 'all') for all
%                           trials in the same plot. Note that 'each' &
%                           scalar plot multiple BDFs on top of each other,
%                           while 'all' generates one figure per BDF
%   (plot_SD)           : [true] plot the data from each target in the
%                           same subplot
%
%
% Note: this code is pretty messy, it'd be good to clean it up, overall the
% last part
%
%

function plot_single_trial_data( single_trial_data, var, params, varargin )


% -------------------------------------------------------------------------
% read inputs
if nargin == 4
    target              = varargin{1};
    plot_SD             = false;
elseif nargin == 5
    target              = varargin{1};
    plot_SD             = varargin{2};
else
    target              = 'each';
    plot_SD             = false;
end



% -------------------------------------------------------------------------
% get some info about what to plot 


% nbr of tasks (they will be plotted on top of each other)
nbr_bdfs                = length(single_trial_data);
% --if single_trial_data is not a cell of tasks but a single task, convert
% it to a cell, so the code doesn't break
if nbr_bdfs == 1
    single_trial_data_new{1} = single_trial_data;
    clear single_trial_data;
    single_trial_data   = single_trial_data_new;
    clear single_trial_data_new;
end


% nbr of vars (neurons, forces, ...)
nbr_vars                = length(params);

% nbr of targets to plot per task
if isnumeric(target)
    
    % if passed the target number
    targets_to_plot     = cell(1,nbr_bdfs);
    for b = 1:nbr_bdfs
        targets_to_plot{b} = target;
    end
    % if numel( unique(cell2mat(targets_to_plot)) ) > 1, error('plotting multiple trials not implemented yet!'); end
    nbr_targets_to_plot = length(target);
    
elseif ischar(target)
    switch target
        % if want to plot each target (in a separate plot)
        case 'each'
            % create a cell with nbr_bdfs fields, each containing a vector
            % with the targets to plot for that BDF
            targets_to_plot = arrayfun( @(y) 1:y, cellfun( @(x) length(x.target)-1, single_trial_data ), ...
                                'UniformOutput', false );
            % get maximum number of targets across all tasks
            nbr_targets_to_plot = length( 1:max(cellfun( @(x) length(x.target)-1, ...
                                    single_trial_data )) );
                                
        % if want to plot all concatenated trials (last target in
        % single_trial_data struct)
        case 'all'
            % create a cell with nbr_bdfs fields, each containing the
            % position of the last target (i.e. the concatenated trials
            % across all targets)
            if ~iscell(single_trial_data)
                targets_to_plot = num2cell( arrayfun( @(x) 1:length(x.target)-1, single_trial_data, ...
                                    'UniformOutput', false ), 1 );
            else
                targets_to_plot = cellfun( @(x) 1:length(x.target)-1, single_trial_data, ...
                                    'UniformOutput', false );
            end
                
            % only plotting one target per task now:
            nbr_targets_to_plot = length(targets_to_plot{1});
    end
end


% ---------------------------------
% what variable to plot
switch var
    case 'neuron'
        var_type{1}     = 'neural_data';
        var_type{2}     = 'fr';
    case 'emg'
        var_type{1}     = 'emg_data';
        var_type{2}     = 'emg';
    case 'pos'
        var_type{1}     = 'pos';
        var_type{2}     = 'data';
    case 'vel'
        var_type{1}     = 'vel';
        var_type{2}     = 'data';
    case 'force'
        var_type{1}     = 'force';
        var_type{2}     = 'data';
    case 'neural PC'
        var_type{1}     = 'neural_data';
        var_type{2}     = 'dim_red';
        var_type{3}     = 'st_scores';
    case 'muscle synergy'
        var_type{1}     = 'emg_data';
        var_type{2}     = 'dim_red';
        var_type{3}     = 'st_scores';
end
        

% ---------------------------------
% define colors & nbr rows/cols

if ischar(target)
    % for the case that we want to plot all the targets one on top of
    % another...
    if strcmp(target,'all')
        colors_plot     = parula(nbr_targets_to_plot);
    end
end
%
if ~exist('colors_plot','var')
    colors_plot         = parula(nbr_vars*nbr_bdfs);
end

% rows and columns in the plot
cols_plot               = ceil(sqrt(nbr_vars));
rows_plot               = ceil(nbr_vars/cols_plot);


% ---------------------------------
% create time vectors
nbr_bins_p_bdf          = cellfun( @(x) size(x.target{1}.neural_data.fr,1), ...
                            single_trial_data );
for b = 1:nbr_bdfs
    t_axis{b}           = 0 : single_trial_data{b}.target{1}.bin_size : ...
                            single_trial_data{b}.target{1}.bin_size*(nbr_bins_p_bdf(b)-1);
end
    


% -------------------------------------------------------------------------
% Plot
    
% get list of targets --can be different across tasks, and whether we want
% to plot all concatenated targets or just one of them
possible_tgts           = cellfun( @(x) 1:length(x.target), single_trial_data, ...
                            'UniformOutput', false );


% see if we want everything in a single figure or one per target

if ischar(target)

    % ---------------------------------------------------------------
    % for the case that we want to plot all the targets one on top of
    % another...
    if strcmp(target,'all')
        for b = 1:nbr_bdfs
            % one figure per BDF
            figure
            
            % independent plot in this figure for each neuron/force/emg...
            for v = 1:nbr_vars

                % retrieve subplot
                subplot(rows_plot,cols_plot,v), hold on

                % now plot mean +/- SD per target in this plot
                for t = 1:nbr_targets_to_plot

                    if ismember(targets_to_plot{b}(t),possible_tgts{b})
                        this_tgt = targets_to_plot{b}(t);

                        % get mean +/-SD of the variable
                        %
                        % check if var name is a 3-by-1 cell (for the neural pcs or
                        % muscle synergies) or a 2-by-1 cell (for the rest)
                        if length(var_type) == 2
                            aux_data    = squeeze(single_trial_data{b}.target{this_tgt}.(var_type{1}).(var_type{2})...
                                (:,params(v),:));
                        else
                            aux_data    = squeeze(single_trial_data{b}.target{this_tgt}.(var_type{1}).(var_type{2}).(var_type{3})...
                                (:,params(v),:));
                        end
                        aux_mean        = mean(aux_data,2);
                        if plot_SD
                            aux_sd      = std(aux_data,0,2);
                        end
                        color_indx      = t;

                        plot( t_axis{b}, aux_mean, 'color',colors_plot(color_indx,:),'linewidth',3);
                        if plot_SD
                            plot( t_axis{b}, aux_mean+aux_sd, ':','color',colors_plot(color_indx,:),'linewidth',3);
                            plot( t_axis{b}, aux_mean-aux_sd, ':','color',colors_plot(color_indx,:),'linewidth',3);
                        end
                        set(gca,'TickDir','out'), set(gca,'FontSize',14);
                        
                        % add labels
                        if v > cols_plot*(rows_plot-1)
                            xlabel('time (s)')
                        end
                        if strcmp(var,'emg')
                            ylabel(single_trial_data{1}.target{1}.emg_data.emg_names{v});
                        else
                            ylabel([var ' ' num2str(v)],'Interpreter','none');
                        end
                    end
                end
            end
        end
    else
        % dirty fix
        aux_each_target     = true; 
    end
end
% -------------------------------------------------------------------
% when we want one figure per target
if isnumeric(target) || exist('aux_each_target','var')
    for t = 1:nbr_targets_to_plot

        % independent plot per target
        % figure('units','normalized','outerposition',[0 0 1 1])
        figure

        % overlay plot per BDF (task)
        for b = 1:nbr_bdfs
            % independent plot in this figure for each neuron/force/emg...
            for v = 1:nbr_vars

                if ismember(targets_to_plot{b}(t),possible_tgts{b})
                    this_tgt        = targets_to_plot{b}(t);

                    % the corresponding subplot
                    subplot(rows_plot,cols_plot,v), hold on

                    % get mean +/-SD of the variable 
                    %
                    % check if var name is a 3-by-1 cell (for the neural pcs or
                    % muscle synergies) or a 2-by-1 cell (for the rest) 
                    if length(var_type) == 2
                        aux_data    = squeeze(single_trial_data{b}.target{this_tgt}.(var_type{1}).(var_type{2})...
                                        (:,params(v),:));
                    else
                        aux_data    = squeeze(single_trial_data{b}.target{this_tgt}.(var_type{1}).(var_type{2}).(var_type{3})...
                                        (:,params(v),:));
                    end
                    aux_mean        = mean(aux_data,2);
                    if plot_SD
                        aux_sd      = std(aux_data,0,2);
                    end
                    color_indx      = v+(b-1)*nbr_vars;
                    if nbr_bdfs == 1, plot( t_axis{b}, aux_data, 'color',[.65 .65 .65]); end
                    plot( t_axis{b}, aux_mean, 'color',colors_plot(color_indx,:),'linewidth',3);
                    if plot_SD
                        plot( t_axis{b}, aux_mean+aux_sd, ':','color',colors_plot(color_indx,:),'linewidth',3);
                        plot( t_axis{b}, aux_mean-aux_sd, ':','color',colors_plot(color_indx,:),'linewidth',3);
                    end
                    set(gca,'TickDir','out'), set(gca,'FontSize',14);
                end

                % add labels
                if v > cols_plot*(rows_plot-1)
                    xlabel('time (s)')
                end
                ylabel([var ' ' num2str(v)],'Interpreter','none');
            end
        end
    end
end
