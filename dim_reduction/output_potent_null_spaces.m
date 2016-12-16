%
% Finf output potent and null spaces of a neural manifold
%
%   function
%
%
%


function opn_spaces = output_potent_null_spaces( single_trial_data, varargin )

% -------------------------------------------------------------------------
% read inputs

if nargin == 2
    params                  = output_potent_null_spaces_defaults( varargin{1} );
else
    params                  = output_potent_null_spaces_defaults();
end


% -------------------------------------------------------------------------
% some checks

% check that the output var is 'emg','pos','vel' or 'force'
output_vars_opts            = {'emg','pos','vel','force','dim_red_emg'};

if sum(strcmpi(output_vars_opts,params.output)) == 0
    error('The output variable for fitting a model can only be emg, pos, vel or force')
end

% see that the neural to output-var delay is a multiple of the bin_size,
% otherwise round it to the smallest integer multiple
if abs( rem(params.neural_to_output_delay,single_trial_data.target{1}.bin_size) ) > 1E-6
    params.neural_to_output_delay = floor(params.neural_to_output_delay/single_trial_data.target{1}.bin_size) * ...
                                    single_trial_data.target{1}.bin_size;
    warning(['neural to ' params.output ' delay is not an exact multiple of the bin size. Set it to: ' ...
        num2str(params.neural_to_output_delay) ' s']);
end



% -------------------------------------------------------------------------
% prepare the data for the regression


% convert input-output lag to nbr of bins -- using 'floor' to avoid
% rounding errors
bins_lag                    = floor(params.neural_to_output_delay/single_trial_data.target{1}.bin_size);

% get target number
switch params.target
    case 'all_conc'
        tgt                 = length(single_trial_data.target);
    otherwise
        tgt                 = params.target;
end


% ----------------------------------------------
% get input data - spikes or dim_red neural data
switch params.input
    case 'spikes'
        % for all (concatenated trials)
        if ~params.trial_averaged
            input_data      = single_trial_data.target{tgt}.neural_data.fr;
        % for the trial-averaged firing rates (PSTH's)
        else
            % create matrix with same format as the single trial data
            bins_p_trial    = length(single_trial_data.target{1}.t);
            input_data      = zeros(bins_p_trial, ...
                                size(single_trial_data.target{tgt}.neural_data.mn,2), ...
                                tgt-1);
            % and fill it with each PSTH
            for t = 1:tgt-1
                input_data(:,:,t) = single_trial_data.target{tgt}.neural_data.mn( ...
                                bins_p_trial*(t-1)+1 : bins_p_trial*t, : );
            end
        end
    case 'dim_red'
         % for all (concatenated trials)
        if ~params.trial_averaged 
            input_data      = single_trial_data.target{tgt}.neural_data.dim_red.st_scores...
                                (:,1:params.dim_neural_manifold,:);
        % for the trial-averaged firing rates (PSTH's)
        else
             % create matrix with same format as the single trial data
            bins_p_trial    = length(single_trial_data.target{1}.t);
            input_data      = zeros(bins_p_trial, ...
                                params.dim_neural_manifold, ...
                                tgt-1);
            % and fill it with each PSTH
            for t = 1:tgt-1
                input_data(:,:,t) = single_trial_data.target{tgt}.neural_data.dim_red.st_scores_mn( ...
                                bins_p_trial*(t-1)+1 : bins_p_trial*t, 1:params.dim_neural_manifold );
            end
        end
end


% ----------------------------------------------
% get output data - emg, pos, vel, force, dim red emg
switch params.output
    case 'emg'
         % for all (concatenated trials)
        if ~params.trial_averaged 
            output_data     = single_trial_data.target{tgt}.emg_data.emg;
        % for the trial-averaged firing rates (PSTH's)
        else
            % create matrix with same format as the single trial data
            bins_p_trial    = length(single_trial_data.target{1}.t);
            output_data     = zeros(bins_p_trial, ...
                                size(single_trial_data.target{tgt}.emg_data.mn,2), ...
                                tgt-1);
            % and fill it with each PSTH
            for t = 1:tgt-1
                output_data(:,:,t) = single_trial_data.target{tgt}.emg_data.mn( ...
                                bins_p_trial*(t-1)+1 : bins_p_trial*t, : );
            end
        end
    case {'pos','vel','force'} % these 2 have the same data structure...
        % for all (concatenated trials)
        if ~params.trial_averaged 
            output_data     = single_trial_data.target{tgt}.pos.data;
        % for the trial-averaged firing rates (PSTH's)
        else
            % create matrix with same format as the single trial data
            bins_p_trial    = length(single_trial_data.target{1}.t);
            output_data     = zeros(bins_p_trial, ...
                                size(single_trial_data.target{tgt}.(params.output).mn,2), ...
                                tgt-1);
            % and fill it with each PSTH
            for t = 1:tgt-1
                output_data(:,:,t) = single_trial_data.target{tgt}.(params.output).mn( ...
                                bins_p_trial*(t-1)+1 : bins_p_trial*t, : );
            end
        end
    case 'dim_red_emg'
        % for all (concatenated trials)
        if ~params.trial_averaged 
            switch single_trial_data.target{end}.emg_data.dim_red.method
                % for PCA, take the nbr of projections defined in
                % params.dim_emg_manifold
                case 'pca'
                    output_data = single_trial_data.target{tgt}.emg_data.dim_red.st_scores...
                                (:,1:params.dim_emg_manifold,:);
                % For NMF, take all that have been defined
                case 'nmf'
                    output_data = single_trial_data.target{tgt}.emg_data.dim_red.st_scores;
            end
        % for the trial-averaged firing rates (PSTH's)
        else
            % create matrix with same format as the single trial data
            bins_p_trial    = length(single_trial_data.target{1}.t);
            switch single_trial_data.target{end}.emg_data.dim_red.method
                % for PCA, take the nbr of projections defined in
                % params.dim_emg_manifold
                case 'pca'
                    output_data  = zeros(bins_p_trial, ...
                                params.dim_emg_manifold, ...
                                tgt-1);
                % For NMF, take all that have been defined                            
                case 'nmf'
                    output_data  = zeros(bins_p_trial, ...
                                size(single_trial_data.target{end}.emg_data.dim_red.w,2), ...
                                tgt-1);
            end
            % and fill it with each PSTH
            for t = 1:tgt-1
                output_data(:,:,t) = single_trial_data.target{tgt}.emg_data.dim_red.st_scores_mn( ...
                                bins_p_trial*(t-1)+1 : bins_p_trial*t, 1:size(output_data,2) );
            end
        end
end



% align them with the desired delay to do linear regression
[X, y]                      = concatenate_to_regress( input_data, output_data, bins_lag );


% check that X has full rank, otherwise give a warning --this would
% indicate that it would be good to cross-validate
if rank(X) < size(X,2)
    warning('matrix X does not have full rank');
end


% -------------------------------------------------------------------------
% Do linear regression to find the neuron to output model
%
% The model will have the form: M = W·N, where
%   M is an m x t matrix with outputs (e.g., EMGs) versus time
%   N is an n x t matrix with inputs (e.g.m spikes) versus time
%   W is an m x n matrix with MISO linear models as rows


% get nbr of outputs
nbr_outputs                 = size(output_data,2);


% Build the cross-validated model W
% -- W_offset is a matrix with the offsets in the model, and stats is a
% matrix with nbr_outpus rows, and 
[W, W_offset, stats]        = build_lm( X, y, true, params.fold_length );


% -------------------------
% compute model predictions
y_hat                       = W*X' + repmat(W_offset,1,size(y,1));
y_hat                       = y_hat';

% -------------------------------------------------------------------------
% Find potent and null spaces by doing SVD of the matrix W
%
% W = U·S·V' where V' defines the task potent space -- note that svd
% returns V not V'

% SVD
[U, S, V]                   = svd( W );

% The output potent spaces is defined by the first m columns of V, where m
% is the number of dimensions of the output
% V_potent                    = V(1:nbr_outputs,:)';
% V_null                      = V(nbr_outputs+1:end,:)';
% V_potent                    = V(:,1:nbr_outputs);
V_potent                    = V(:,1:nbr_outputs);
V_null                      = V(:,nbr_outputs+1:end);


% -------------------------------------------------------------------------
% Return vars

opn_spaces.svdec.U          = U;
opn_spaces.svdec.S          = S;
opn_spaces.svdec.V          = V;
opn_spaces.W                = W;
opn_spaces.V_potent         = V_potent;
opn_spaces.V_null           = V_null;
opn_spaces.stats_W          = stats;
opn_spaces.lm_fit.X         = X;
opn_spaces.lm_fit.y         = y;
opn_spaces.lm_fit.y_hat     = y_hat;
opn_spaces.lm_fit.W_offset  = W_offset;