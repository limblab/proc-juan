%
%
%

function opn_results = batch_compare_potent_null_spaces( varargin ) 


% check that the parallel pool is running, otherwise start it
gcp;

% -------------------------------------------------------------------------
% read inputs

if nargin >= 1
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end
if nargin >= 2
    params              = batch_compare_potent_null_spaces_defaults( varargin{2:end} );
else
    params              = batch_compare_potent_null_spaces_defaults();
end


% -------------------------------------------------------------------------
% load data

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end


% -------------------------------------------------------------------------
% get info about monkeys, tasks, and who did what

meta_info               = batch_get_monkey_task_data( datasets );


% -------------------------------------------------------------------------
% compute output potent and null spaces for each task in each session, and
% compare them

data                    = cell(1,length(meta_info.tasks_per_session));

for i = 1:meta_info.nbr_monkeys

    for ii = 1:meta_info.nbr_sessions_per_monkey(i)
       
        % what dataset are we looking at?
        dtst            = meta_info.sessions_per_monkey{i}(ii);

        disp(['comparing output potent / null spaces dataset #' num2str(dtst)]);

        
        % --------------------------------------
        % compute output potent / null spaces
        
        % choose the time window
        params_orig     = params;
        params.time_win = params.time_win(dtst,:);
        
        % and do the analysis
        opn_spaces      = output_potent_null_spaces_all_manifolds( datasets{dtst}.stdata, ...
                            params );
                        
        % revert params
        params          = params_orig;
        
        
        % --------------------------------------
        % compare output potent spaces
        opn_space_comps = compare_output_potent_null_spaces( opn_spaces );
        
        
        % --------------------------------------
        % store the results
        data{dtst}.opn_spaces       = opn_spaces;
        data{dtst}.opn_space_comps  = opn_space_comps;
    end
end




% -------------------------------------------------------------------------
% return data struct

% generate return struct
opn_results.data        = data; 





% -------------------------------------------------------------------------
% plots

% some quick & dirty plotting

% Plot the R2 of the model fits
colors                  = distinguishable_colors(length(datasets));
ctr                     = 1;
figure, hold on
for i = 1:length(datasets)
    for ii = 1:length(data{i}.opn_spaces)
        plot(ctr,data{i}.opn_spaces(ii).stats_W(:,1),'.','color',colors(i,:),'markersize',12) 
        ctr             = ctr + 1;
    end
end
set(gca,'Tickdir','out'),set(gca,'FontSize',14)
ylim([0 1]),xlim([0, ctr]),xlabel('task #'),ylabel('R^2')


% generate random principal angle distribution, to compare
plane_dim               = unique(cellfun(@(x) size(x.opn_spaces(1).V_potent,2), opn_results.data ));
plane_dim(plane_dim>params.dim_neural_manifold) = [];
[~,angle_non_orth]      = empirical_principal_angle_distribution( params.dim_neural_manifold, plane_dim, ...
                            100000, params.P_thr, false );

all_theta               = cellfun(@(x) x.opn_space_comps.princ_ang_V_potent, data,...
                            'UniformOutput', false );

% plot empirical angle distributions for each task
for i = 1:length(datasets)
   figure, hold on
   plot(rad2deg(all_theta{i}'),'color',colors(i,:),'linewidth',2) 
   indx_this            = find( plane_dim == size(opn_results.data{1}.opn_spaces(1).V_potent,2) );
   plot(angle_non_orth(1:size(opn_results.data{1}.opn_spaces(1).V_potent,2),indx_this),...
        'color',[.6 .6 .6],'linewidth',3,'linestyle',':')
   set(gca,'Tickdir','out'),set(gca,'FontSize',14)
   xlabel('dimension'),ylabel('principal angle (deg)')
   xlim([0, size(opn_results.data{1}.opn_spaces(1).V_potent,2)+1]),ylim([0 90])
end


