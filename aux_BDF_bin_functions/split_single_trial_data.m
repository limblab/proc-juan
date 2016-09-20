%
% This function splits a single_trial_data struct in two, for doing various
% control analyses
%
%   function single_trial_data = split_single_trial_data( single_trial_data, varargin )
%
%
% Inputs (opt)          : [default]
%   single_trial_data   : a single single_trial_data struct
%   (nbr_chunks)        : [2] nbr of single_trial_data structs into which
%                           single_trial_data will be split
%   (random_order_yn)   : [false] randomize the order of the trials
%
% Outputs:
%   single_trial_data   : cell array w the split single_trial_data structs
%
%


function single_trial_data = split_single_trial_data( single_trial_data, varargin )



% -------------------------------------------------------------------------
% read inputs
if nargin >= 2
    nbr_chunks          = varargin{1};
else
    nbr_chunks          = 2;
end

if nargin == 3
    random_order_yn     = varargin{2};
else
    random_order_yn     = false;
end


% -------------------------------------------------------------------------
% get some meta info

% get nbr of targets
nbr_targets             = length(single_trial_data.target)-1;

% get nbr of trials_p_target
nbr_trials_p_target     = cellfun( @(x) size(x.bin_indx_p_trial,2), single_trial_data.target );
nbr_trials_p_target(end) = [];

% compute nbr of trials to keep for each target
nbr_trials_to_keep_p_target = floor(nbr_trials_p_target/nbr_chunks);


% offset to know where to start counting for each trial
trial_offset            = zeros(1,nbr_targets);
aux_to                  = repmat(nbr_trials_to_keep_p_target,nbr_chunks-1,1);
trial_offset            = [trial_offset; cumsum(aux_to,1)];
clear aux_to;


% create a cell array with the trials to keep for each target and "chunk"
for t = 1:nbr_targets
    % if want to keep consecutive trials
    if ~random_order_yn
        % get start and end trial indexes for the n trials that each chunk
        % will comprise of --note that they can be different for each
        % target, as the number of trials does not necessarily have to be
        % the same
        for c = 1:nbr_chunks
            start_this      = 1 + trial_offset(c,t);
            end_this        = nbr_trials_to_keep_p_target(t) + trial_offset(c,t);
            indxs_this{c,t} = start_this:end_this;
        end
    % or randomly ordered_trials
    else
        % randomly order the trials and then take the n trials that will be
        % each chunk
        aux_rand_indexes    = randperm(nbr_trials_p_target(t));
        for c = 1:nbr_chunks
            start_this      = 1 + trial_offset(c,t);
            end_this        = nbr_trials_to_keep_p_target(t) + trial_offset(c,t);
            indxs_this{c,t} = aux_rand_indexes(start_this:end_this);
        end
    end
end


% create a cell array with the trials to keep for each target and "chunk"
for t = 1:nbr_targets

    % if want to keep consecutive trials just keep the number that
    % correspond to each chunk
    if ~random_order_yn
        aux_rand_indexes = 1:nbr_trials_p_target(t);
    % if we want to randomly choose the trials for each chunk
    else
        aux_rand_indexes = randperm(nbr_trials_p_target(t));
    end

    % store them here
    for c = 1:nbr_chunks
        start_this      = 1 + trial_offset(c,t);
        end_this        = nbr_trials_to_keep_p_target(t) + trial_offset(c,t);
        indxs_this{c,t} = aux_rand_indexes(start_this:end_this);
    end
end


% see what optional fields we have
smoothed_fr_yn          = isfield(single_trial_data.target{1}.neural_data,'smoothed_fr');
dim_red_neural_yn       = isfield(single_trial_data.target{1}.neural_data,'dim_red');
dim_red_emg_yn          = isfield(single_trial_data.target{1}.emg_data,'dim_red');
pos_yn                  = isfield(single_trial_data.target{1},'pos');
vel_yn                  = isfield(single_trial_data.target{1},'vel');
force_yn                = isfield(single_trial_data.target{1},'force');

% cell array to go through pos, vel & force in a loop, as the structures
% are identical
other_vars              = {'pos','vel','force'};
other_vars_yn           = [pos_yn, vel_yn, force_yn];

% cell array to go through some meta info fields that will just be copied
% from the original struct to each new struct
meta_fields             = {'bin_size','avg_method','w_i','w_f','t'};


% -------------------------------------------------------------------------
% and split the data

stdata                  = cell(1,nbr_chunks);

% do for the individual targets --but not for the {end} one yet
for c = 1:nbr_chunks
    
    for t = 1:nbr_targets
    
        
        % -----------------------------------------------------------------
        % store the data
        
        % -----------------------
        % 1) neural data vars
        
        % a) FRs
        stdata{c}.target{t}.neural_data.fr              = single_trial_data.target{t}.neural_data.fr(:,:,indxs_this{c,t});
        
        if smoothed_fr_yn
            stdata{c}.target{t}.neural_data.smoothed_fr = single_trial_data.target{t}.neural_data.smoothed_fr(:,:,indxs_this{c,t});
        end
        
        % b) summary stats
        stdata{c}.target{t}.neural_data.mn              = mean(stdata{c}.target{t}.neural_data.fr,3);
        stdata{c}.target{t}.neural_data.sd              = std(stdata{c}.target{t}.neural_data.fr,0,3);
        
        if smoothed_fr_yn
            stdata{c}.target{t}.neural_data.smoothed_fr_mn = mean(stdata{c}.target{t}.neural_data.smoothed_fr,3);
            stdata{c}.target{t}.neural_data.smoothed_fr_sd = std(stdata{c}.target{t}.neural_data.smoothed_fr,0,3);
        end
        
        % c) channel nbrs
        stdata{c}.target{t}.neural_data.neural_chs      = single_trial_data.target{t}.neural_data.neural_chs;
        
        
        % d) dim_red FRs
        if dim_red_neural_yn
            % single trial data
            stdata{c}.target{t}.neural_data.dim_red.st_scores = ...
                single_trial_data.target{t}.neural_data.dim_red.st_scores(:,:,indxs_this{c,t});
            % summary stats
            stdata{c}.target{t}.neural_data.dim_red.st_scores_mn = mean(stdata{c}.target{t}.neural_data.dim_red.st_scores,3);
            stdata{c}.target{t}.neural_data.dim_red.st_scores_sd = std(stdata{c}.target{t}.neural_data.dim_red.st_scores,0,3);
        end
        
        
        % -----------------------
        % emg data vars
        
        % a) EMGs
        stdata{c}.target{t}.emg_data.emg                = single_trial_data.target{t}.emg_data.emg(:,:,indxs_this{c,t});
        
        % b) summary stats
        stdata{c}.target{t}.emg_data.mn                 = mean(stdata{c}.target{t}.emg_data.emg,3);
        stdata{c}.target{t}.emg_data.sd                 = std(stdata{c}.target{t}.emg_data.emg,0,3);
        
        % c) emg channels
        stdata{c}.target{t}.emg_data.emg_names          = single_trial_data.target{t}.emg_data.emg_names;
        
        
        % d) dim_red EMGs
        if dim_red_emg_yn
            % single trial muscle synergies
            stdata{c}.target{t}.emg_data.dim_red.st_scores = single_trial_data.target{t}.emg_data.dim_red.st_scores(:,:,indxs_this{c,t});
            % summary stats
            stdata{c}.target{t}.emg_data.dim_red.st_scores_mn = mean(stdata{c}.target{t}.emg_data.dim_red.st_scores,3);
            stdata{c}.target{t}.emg_data.dim_red.st_scores_sd = std(stdata{c}.target{t}.emg_data.dim_red.st_scores,0,3);
        end
        
        
        % -----------------------
        % pos, vel & force data
        for v = 1:length(other_vars)
           
            if other_vars_yn(v)
                % single trial pos/vel/force
                stdata{c}.target{t}.(other_vars{v}).data = single_trial_data.target{t}.(other_vars{v}).data(:,:,indxs_this{c,t});
                % summary stats
                stdata{c}.target{t}.(other_vars{v}).mn  = mean( stdata{c}.target{t}.(other_vars{v}).data, 3);
                stdata{c}.target{t}.(other_vars{v}).sd  = std( stdata{c}.target{t}.(other_vars{v}).data, 0, 3);
            end
        end
        
        
        % -----------------------
        % some other stuff
        
        % bin indx
        stdata{c}.target{t}.bin_indx_p_trial            = single_trial_data.target{t}.bin_indx_p_trial(:,indxs_this{c,t});
        
        % bin_size, method for getting the trials (cutting vs. time
        % warping), words for finding the trials, and time vector
        for v = 1:length(meta_fields)
            stdata{c}.target{t}.(meta_fields{v})        = single_trial_data.target{t}.(meta_fields{v});
        end
    end
end



% -------------------------------------------------------------------------
% add last target, which is all the concatenated targets
for c = 1:length(stdata)
    stdata{c}           = add_conc_target_single_trial_data( stdata{c} );
end

% add a some info about dim_red_emg, dim_red_FR
for c = 1:length(stdata)

    if dim_red_neural_yn
        stdata{c}.target{end}.neural_data.dim_red.method    = single_trial_data.target{end}.neural_data.dim_red.method;
        stdata{c}.target{end}.neural_data.dim_red.w         = single_trial_data.target{end}.neural_data.dim_red.w;
        stdata{c}.target{end}.neural_data.dim_red.eigen     = single_trial_data.target{end}.neural_data.dim_red.eigen;
        stdata{c}.target{end}.neural_data.dim_red.chs       = single_trial_data.target{end}.neural_data.dim_red.chs;
    end
    
    if dim_red_emg_yn
        stdata{c}.target{end}.emg_data.dim_red.method       = single_trial_data.target{end}.emg_data.dim_red.method;
        stdata{c}.target{end}.emg_data.dim_red.w            = single_trial_data.target{end}.emg_data.dim_red.w;
        stdata{c}.target{end}.emg_data.dim_red.chs          = single_trial_data.target{end}.emg_data.dim_red.chs;
        if isfield(single_trial_data.target{end}.emg_data,'eigen')
            stdata{c}.target{end}.emg_data.dim_red.eigen    = single_trial_data.target{end}.emg_data.dim_red.eigen;
        end
    end
end


% -------------------------------------------------------------------------
% add concatenated fields (conc_t, conc_fr, conc_emg, ...) to each target
% in the struct

for c = 1:length(stdata)

    stdata{c}           = concatenate_single_trials( stdata{c} );
end


% -------------------------------------------------------------------------
% return var
clear single_trial_data
single_trial_data       = stdata;