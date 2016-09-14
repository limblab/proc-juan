%
% concatenate single trial neural, emg, pos, vel, force data from a
% single_trial_data struct and add them as a "concatenated" time x chs
% array to the respective field in the single_trial_data struct 
%
%   function single_trial_data = concatenate_single_trials( single_trial_data )
%
%


function single_trial_data = concatenate_single_trials( single_trial_data )


% do for all the targets, including the one that is all the concatenated
% trials
for t = 1:length(single_trial_data.target)

    
    % reshape the 3-D time-by-channel-by-trial matrices to 2-D time-by-channel 
    % matrices, with all the trials concatenated
    
    % 1) the FRs
    aux                 = single_trial_data.target{t}.neural_data.fr;
    conc                = zeros(size(aux,1)*size(aux,3),size(aux,2));    
    for d = 1:size(aux,2)
        sq_conc         = reshape( squeeze(aux(:,d,:)), [], 1 );
        conc(:,d)       = sq_conc;
    end
    
    single_trial_data.target{t}.neural_data.conc_fr = conc;
    
    
    % 2) the smoothed FRs --if present
    if isfield(single_trial_data.target{t}.neural_data,'smoothed_fr')
        aux             = single_trial_data.target{t}.neural_data.smoothed_fr;
        conc            = zeros(size(aux,1)*size(aux,3),size(aux,2));    
        for d = 1:size(aux,2)
            sq_conc     = reshape( squeeze(aux(:,d,:)), [], 1 );
            conc(:,d)   = sq_conc;
        end
        
        single_trial_data.target{t}.neural_data.conc_smoothed_fr = conc;
    end

    
    % 3) the dimensionally-reduced FRs --if present
    if isfield(single_trial_data.target{t}.neural_data,'dim_red')
        aux             = single_trial_data.target{t}.neural_data.dim_red.st_scores;
        conc            = zeros(size(aux,1)*size(aux,3),size(aux,2));    
        for d = 1:size(aux,2)
            sq_conc     = reshape( squeeze(aux(:,d,:)), [], 1 );
            conc(:,d)   = sq_conc;
        end
        
        single_trial_data.target{t}.neural_data.dim_red.scores = conc;
    end
    
    
    % 4) the EMGs
    aux                 = single_trial_data.target{t}.emg_data.emg;
    conc                = zeros(size(aux,1)*size(aux,3),size(aux,2));    
    for d = 1:size(aux,2)
        sq_conc         = reshape( squeeze(aux(:,d,:)), [], 1 );
        conc(:,d)       = sq_conc;
    end
    
    single_trial_data.target{t}.emg_data.conc_emg = conc;
    
    
    % 5) the dimensionally-reduced EMGs
    if isfield(single_trial_data.target{t}.emg_data,'dim_red')
        
        aux             = single_trial_data.target{t}.emg_data.dim_red.scores;
        conc            = zeros(size(aux,1)*size(aux,3),size(aux,2));    
        for d = 1:size(aux,2)
            sq_conc     = reshape( squeeze(aux(:,d,:)), [], 1 );
            conc(:,d)   = sq_conc;
        end
        
        single_trial_data.target{t}.emg_data.dim_red.scores = conc;
    end
    
    
    % 6) the pos and vel
    if isfield(single_trial_data.target{t},'pos')
        
        % position
        aux             = single_trial_data.target{t}.pos.data;
        conc            = zeros(size(aux,1)*size(aux,3),size(aux,2));
        for d = 1:size(aux,2)
            sq_conc     = reshape( squeeze(aux(:,d,:)), [], 1 );
            conc(:,d)   = sq_conc;
        end
        
        single_trial_data.target{t}.pos.conc_data = conc;
        
        % velocity
        aux             = single_trial_data.target{t}.vel.data;
        conc            = zeros(size(aux,1)*size(aux,3),size(aux,2));
        for d = 1:size(aux,2)
            sq_conc     = reshape( squeeze(aux(:,d,:)), [], 1 );
            conc(:,d)   = sq_conc;
        end
        
        single_trial_data.target{t}.vel.conc_data = conc;
    end
    

    % add concatenated (continuous) time (arbitrarily starting at 0)
    single_trial_data.target{t}.conc_t = 0:single_trial_data.target{t}.bin_size:...
        (size(single_trial_data.target{t}.neural_data.conc_fr,1)-1)*...
        single_trial_data.target{t}.bin_size;
    % add concatenated bin indexes
    single_trial_data.target{t}.bin_indx = reshape(single_trial_data.target{t}.bin_indx_p_trial,1,[]);
    
%     % overwrite the time vector associated to the "neural scores" keeping
%     % its integrity
%     if isfield(single_trial_data.target{t}.neural_data,'dim_red')
%     
%         aux_indx_trans      = find(diff(single_trial_data.target{t}.neural_data.dim_red.t) > ...
%                                 single_trial_data.target{t}.bin_size*1.1);
%         aux_indx_trans      = [1; aux_indx_trans+1];
%         indx_dim_red_neural_t = zeros( size(single_trial_data.target{t}.neural_data.dim_red.st_scores,3),...
%                                 size(single_trial_data.target{t}.neural_data.dim_red.st_scores,1) );
%         for ii = 1:length(aux_indx_trans)
%             indx_dim_red_neural_t(ii,:) = aux_indx_trans(ii) : aux_indx_trans(ii) + ...
%                                 size(single_trial_data.target{t}.neural_data.dim_red.st_scores,1) - 1;
%         end
%         indx_dim_red_neural_t = sort(reshape(indx_dim_red_neural_t,1,[]));

%         single_trial_data.target{t}.neural_data.dim_red.t = ...
%             single_trial_data.target{t}.neural_data.dim_red.t(indx_dim_red_neural_t,1);
    end
    
end
