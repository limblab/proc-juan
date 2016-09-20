%
% Add an additional 'target' field to a single_trial_data struct, that has
% the ata from all the trials together
%
%   single_trial_data = add_conc_target_single_trial_data( single_trial_data )
%
%

function single_trial_data = add_conc_target_single_trial_data( single_trial_data )


% get nbr of targets
nbr_targets                 = length(single_trial_data.target);

% pointer for the new target
ptr                         = nbr_targets + 1;


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
% Do

% ------------------------------------------
% 1) copy the data from the first target to the new target struct

% --------------
% A) Neural data

% FRs + mean, SD
single_trial_data.target{ptr}.neural_data.fr        = single_trial_data.target{1}.neural_data.fr;

single_trial_data.target{ptr}.neural_data.mn        = single_trial_data.target{1}.neural_data.mn;
single_trial_data.target{ptr}.neural_data.sd        = single_trial_data.target{1}.neural_data.sd;

% smoothed FRs + mean, SD
if smoothed_fr_yn
    single_trial_data.target{ptr}.neural_data.smoothed_fr = single_trial_data.target{1}.neural_data.smoothed_fr;
    
    single_trial_data.target{ptr}.neural_data.smoothed_fr_mn = single_trial_data.target{1}.neural_data.smoothed_fr_mn;
    single_trial_data.target{ptr}.neural_data.smoothed_fr_sd = single_trial_data.target{1}.neural_data.smoothed_fr_sd;
end

% neural chs
single_trial_data.target{ptr}.neural_data.neural_chs = single_trial_data.target{1}.neural_data.neural_chs;


% dim_red neural data
if dim_red_neural_yn
    single_trial_data.target{ptr}.neural_data.dim_red.st_scores = single_trial_data.target{1}.neural_data.dim_red.st_scores;
    
    single_trial_data.target{ptr}.neural_data.dim_red.st_scores_mn = single_trial_data.target{1}.neural_data.dim_red.st_scores_mn;
    single_trial_data.target{ptr}.neural_data.dim_red.st_scores_sd = single_trial_data.target{1}.neural_data.dim_red.st_scores_sd;
end


% --------------
% B) EMG data

% EMG data + mean, SD
single_trial_data.target{ptr}.emg_data.emg          = single_trial_data.target{1}.emg_data.emg;

single_trial_data.target{ptr}.emg_data.mn           = single_trial_data.target{1}.emg_data.mn;
single_trial_data.target{ptr}.emg_data.sd           = single_trial_data.target{1}.emg_data.sd;

% emg nagmes
single_trial_data.target{ptr}.emg_data.emg_names    = single_trial_data.target{1}.emg_data.emg_names;


% dim red EMG data
if dim_red_emg_yn
    
    single_trial_data.target{ptr}.emg_data.dim_red.st_scores = single_trial_data.target{1}.emg_data.dim_red.st_scores;
    
    single_trial_data.target{ptr}.emg_data.dim_red.st_scores_mn = single_trial_data.target{1}.emg_data.dim_red.st_scores_mn;
    single_trial_data.target{ptr}.emg_data.dim_red.st_scores_sd = single_trial_data.target{1}.emg_data.dim_red.st_scores_sd;
end


% --------------
% C) POS/VEL/FORCE data

for v = 1:length(other_vars)
    
    if other_vars_yn(v)
   
        single_trial_data.target{ptr}.(other_vars{v}).data  = single_trial_data.target{1}.(other_vars{v}).data;
        
        single_trial_data.target{ptr}.(other_vars{v}).mn    = single_trial_data.target{1}.(other_vars{v}).mn;
        single_trial_data.target{ptr}.(other_vars{v}).sd    = single_trial_data.target{1}.(other_vars{v}).sd;
    end
end


% --------------
% D) META STUFF + TIME VECTOR (bin_size, avg_method, words, ...)

for v = 1:length(meta_fields)
   
    single_trial_data.target{ptr}.(meta_fields{v})          = single_trial_data.target{1}.(meta_fields{v});
end

single_trial_data.target{ptr}.bin_indx_p_trial              = single_trial_data.target{1}.bin_indx_p_trial;



% ------------------------------------------
% 2) concatenate the data from other trials


for t = 2:nbr_targets
    
    % --------------
    % A) Neural data
    
    % FRs + mean, SD
    single_trial_data.target{ptr}.neural_data.fr        = cat(3,single_trial_data.target{ptr}.neural_data.fr, ...
                                                            single_trial_data.target{t}.neural_data.fr);

    single_trial_data.target{ptr}.neural_data.mn        = cat(1,single_trial_data.target{ptr}.neural_data.mn, ...
                                                            single_trial_data.target{t}.neural_data.mn);
    single_trial_data.target{ptr}.neural_data.sd        = cat(1,single_trial_data.target{ptr}.neural_data.sd, ...
                                                            single_trial_data.target{t}.neural_data.sd);

    % smoothed FRs + mean, SD
    if smoothed_fr_yn
        
        single_trial_data.target{ptr}.neural_data.smoothed_fr = cat( 3, single_trial_data.target{ptr}.neural_data.smoothed_fr, ...
                                                            single_trial_data.target{t}.neural_data.smoothed_fr );

        single_trial_data.target{ptr}.neural_data.smoothed_fr_mn    = cat( 1, single_trial_data.target{ptr}.neural_data.smoothed_fr_mn, ...
                                                                single_trial_data.target{t}.neural_data.smoothed_fr_mn );
        single_trial_data.target{ptr}.neural_data.smoothed_fr_sd    = cat( 1, single_trial_data.target{ptr}.neural_data.smoothed_fr_sd, ...
                                                                single_trial_data.target{t}.neural_data.smoothed_fr_sd );
    end
    
    
    % dim_red neural data
    if dim_red_neural_yn
       
        single_trial_data.target{ptr}.neural_data.dim_red.st_scores = cat( 3, single_trial_data.target{ptr}.neural_data.dim_red.st_scores, ...
                                                                    single_trial_data.target{t}.neural_data.dim_red.st_scores );
                                                                
        single_trial_data.target{ptr}.neural_data.dim_red.st_scores_mn = cat( 1, single_trial_data.target{ptr}.neural_data.dim_red.st_scores_mn, ...
                                                                    single_trial_data.target{t}.neural_data.dim_red.st_scores_mn );
        single_trial_data.target{ptr}.neural_data.dim_red.st_scores_sd = cat( 1, single_trial_data.target{ptr}.neural_data.dim_red.st_scores_sd, ...
                                                                    single_trial_data.target{t}.neural_data.dim_red.st_scores_sd );
    end
    
    
    % --------------
    % B) EMG DATA
    
    % EMG data + mean, SD
    single_trial_data.target{ptr}.emg_data.emg          = cat( 3, single_trial_data.target{ptr}.emg_data.emg, ...
                                                            single_trial_data.target{t}.emg_data.emg );
                                                        
    single_trial_data.target{ptr}.emg_data.mn           = cat( 1, single_trial_data.target{ptr}.emg_data.mn, ...
                                                            single_trial_data.target{t}.emg_data.mn );
    single_trial_data.target{ptr}.emg_data.sd           = cat( 1, single_trial_data.target{ptr}.emg_data.sd, ...
                                                            single_trial_data.target{t}.emg_data.sd );
                                                        
    % dim_red EMG data
    if dim_red_emg_yn
       
        single_trial_data.target{ptr}.emg_data.dim_red.st_scores    = cat( 3, single_trial_data.target{ptr}.emg_data.dim_red.st_scores, ...
                                                            single_trial_data.target{t}.emg_data.dim_red.st_scores );
        single_trial_data.target{ptr}.emg_data.dim_red.st_scores_mn = cat( 1, single_trial_data.target{ptr}.emg_data.dim_red.st_scores_mn, ...
                                                            single_trial_data.target{t}.emg_data.dim_red.st_scores_mn );
        single_trial_data.target{ptr}.emg_data.dim_red.st_scores_sd = cat( 1, single_trial_data.target{ptr}.emg_data.dim_red.st_scores_sd, ...
                                                            single_trial_data.target{ptr}.emg_data.dim_red.st_scores_sd );
    end
    
    
    % --------------
    % C) POS/VEL/FORCE data

    for v = 1:length(other_vars)
    
        if other_vars_yn(v)
            single_trial_data.target{ptr}.(other_vars{v}).data      = cat( 3, single_trial_data.target{ptr}.(other_vars{v}).data, ...
                                                            single_trial_data.target{t}.(other_vars{v}).data );
            single_trial_data.target{ptr}.(other_vars{v}).mn        = cat( 1, single_trial_data.target{ptr}.(other_vars{v}).mn, ...
                                                            single_trial_data.target{ptr}.(other_vars{v}).mn );
            single_trial_data.target{ptr}.(other_vars{v}).sd        = cat( 1, single_trial_data.target{ptr}.(other_vars{v}).sd, ...
                                                            single_trial_data.target{ptr}.(other_vars{v}).sd );

        end
    end
    
    % --------------
    % D) BIN INDX PER TRIAL
    single_trial_data.target{ptr}.bin_indx_p_trial      = cat( 2, single_trial_data.target{ptr}.bin_indx_p_trial, ...
                                                            single_trial_data.target{t}.bin_indx_p_trial ); 
end



