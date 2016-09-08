%
%
%

function batch_compare_potent_null_spaces( varargin ) 


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
        comp_opn_spaces = compare_output_potent_null_spaces( opn_spaces );
        
        
        % --------------------------------------
        % store the results
        data{dtst}.opn_spaces       = opn_spaces;
        data{dtst}.comp_opn_spaces  = comp_opn_spaces;
    end
end




% -------------------------------------------------------------------------
% plots



figure, hold on
