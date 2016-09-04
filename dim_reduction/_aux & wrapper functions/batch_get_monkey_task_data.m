%
% Retrieve information about monkeys in the dataset, the tasks they
% performed, the number of sessions, etc
% 

function meta_info = batch_get_monkey_task_data( varargin )


% -------------------------------------------------------------------------
% read input params

if nargin
    if iscell(varargin{1})
        datasets        = varargin{1};
    else
        path            = varargin{1};
    end
else
    path                = pwd;
end
clear varargin;


% -------------------------------------------------------------------------
% load data if not passed

if ~exist('datasets','var')
    load([path filesep 'all_manifold_datasets.mat'])
end


% -------------------------------------------------------------------------
% Retrieve information about monkeys and tasks

% what monkeys do I have?
meta.monkeys            = unique(cellfun( @(x) x.monkey, datasets, 'UniformOutput', false ));
meta.nbr_monkeys        = numel(meta.monkeys);

% what tasks did they perform in each session?
meta.tasks_per_session  = cellfun( @(x) x.labels, datasets, 'UniformOutput', false );

% what task did each monkey perform?
meta.tasks_per_monkey   = cell(1,meta.nbr_monkeys);
meta.sessions_per_monkey = cell(1,meta.nbr_monkeys);

for i = 1:meta.nbr_monkeys
    datasets_this_monk  = find( cellfun( @(x) strncmp(x.monkey,meta.monkeys{i},...
                            meta.nbr_monkeys), datasets ) );
    meta.tasks_per_monkey{i} = meta.tasks_per_session{datasets_this_monk(1)};
    for ii = 2:length(datasets_this_monk)
        for iii = 1:size(meta.tasks_per_session{datasets_this_monk(ii)},2)
            if sum( strncmp(meta.tasks_per_monkey{i},...
                    meta.tasks_per_session{datasets_this_monk(ii)}(iii),...
                    length(meta.tasks_per_session{datasets_this_monk(ii)}(iii))) ) == 0
                
                meta.tasks_per_monkey{i}(numel(meta.tasks_per_monkey{i})) = ...
                    meta.tasks_per_session{datasets_this_monk(ii)}(iii);
            end
        end
    end
    meta.sessions_per_monkey{i} = datasets_this_monk;
end
meta.nbr_sessions_per_monkey = cell2mat( cellfun( @(x) length(x), meta.sessions_per_monkey, ...
                                'UniformOutput', false ) );

% -------------------------------------------------------------------------
% Make an ordered list of pairs of tasks within each session, session
% number and monkey. All these will be lists of length equal to the number
% of pairs of tasks --this is useful because there are a lot of pairwise
% comparisons in the analyses

nbr_tasks_pairs                 = sum(cellfun(@(x) nchoosek(length(x),2), ...
                                    meta.tasks_per_session ));

meta.task_pairs.monkey          = cell(1,nbr_tasks_pairs);
meta.task_pairs.session         = zeros(1,nbr_tasks_pairs);
meta.task_pairs.task_pair       = cell(1,nbr_tasks_pairs);
ctr                             = 1;

% store the monkey, pair of tasks and session for each pair of tasks
for m = 1:meta.nbr_monkeys
    for s = 1:length(meta.sessions_per_monkey{m})
        
        % find pairs of tasks in this session
        pairs_this_session      = nchoosek(...
            1:numel(meta.tasks_per_session{meta.sessions_per_monkey{m}(s)}),2);
        nbr_pairs_this_session  = size(pairs_this_session,1);
        
        for p = 1:nbr_pairs_this_session
            % save monkey, session number, and current pair of tasks
            meta.task_pairs.monkey{ctr} = meta.monkeys{m};
            meta.task_pairs.session(ctr) = meta.sessions_per_monkey{m}(s);
            
            meta.task_pairs.task_pair{ctr} = meta.tasks_per_session{meta.sessions_per_monkey{m}(s)}(pairs_this_session(p,:)); 
                        
            % sort tasks alphabetically
            meta.task_pairs.task_pair{ctr} = sort(meta.task_pairs.task_pair{ctr});
            
            ctr                 = ctr + 1;
        end
    end
end


% add unique pairs of tasks to the struct
% ~this will find all the pairs of tasks performed by each monkey
meta.task_pairs.unique_pairs    = [];
for m = 1:numel(meta.monkeys)
    indx_pairs_this_monk        = nchoosek(1:size(meta.tasks_per_monkey{m},2),2);
    for p = 1:size(indx_pairs_this_monk,1)
        pairs_this_monk{p}      = meta.tasks_per_monkey{m}(indx_pairs_this_monk(p,:));
    end
    meta.task_pairs.unique_pairs = [meta.task_pairs.unique_pairs, pairs_this_monk];
    pairs_this_monk             = [];
end
%~here we get rid of duplicate pairs
for i = 1:size(meta.task_pairs.unique_pairs,2)
    for ii = i+1:size(meta.task_pairs.unique_pairs,2)
        str_cmp                 = strcmpi( meta.task_pairs.unique_pairs{i}, ...
                                        meta.task_pairs.unique_pairs{ii} );
        % see if both tasks in the pair are the same
        if sum(str_cmp) == 2 
           meta.task_pairs.unique_pairs{ii} = [];
        end
    end
end
% remove the empty 
meta.task_pairs.unique_pairs    = meta.task_pairs.unique_pairs(...
                                        ~cellfun('isempty',meta.task_pairs.unique_pairs));
% and sort each pair alphabetically, to simplify searchs
meta.task_pairs.unique_pairs = cellfun(@(x) sort(x),meta.task_pairs.unique_pairs,...
                                'UniformOutput',false);


% add a number to identify each task pair. This will give an array that
% will make plotting and sorting things easier
meta.task_pairs.task_pair_nbr   = zeros(1,nbr_tasks_pairs);
for i = 1:length(meta.task_pairs.task_pair)
    for ii = 1:length(meta.task_pairs.unique_pairs)
        if sum( strcmpi(meta.task_pairs.task_pair{i},meta.task_pairs.unique_pairs{ii}) ) == 2
            meta.task_pairs.task_pair_nbr(i) = ii;
        end
    end
end

                            
% -------------------------------------------------------------------------
% Some summary stuff

% what tasks do we have? and how many are they?
meta.tasks              = [];
for i = 1:meta.nbr_monkeys
    meta.tasks          = [meta.tasks meta.tasks_per_monkey{i}];
end
meta.tasks              = unique(meta.tasks);
meta.nbr_tasks          = length(meta.tasks);


% return variable
meta_info               = meta;