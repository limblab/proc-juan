%
% Compare output potent and null spaces of input-output (e.g., neurons to
% emg or to force) linear models (with one single bin). The function
% compares all possible pairs of models.
%
%   function comp_opn_spaces = compare_output_potent_null_spaces( opn_spaces )
%
%
% Input:
%   opn_spaces          : struct generated with output_potent_null_spaces,
%                           with each field being a different task
%
% Output:
%   comp_opn_spaecs     : struct with comparison results
%
%

function comp_opn_spaces = compare_output_potent_null_spaces( opn_spaces ) 


% -------------------------------------------------------------------------
% get nbr of tasks
nbr_bdfs                = length(opn_spaces);
% and all possible combinations
comb_bdfs               = nchoosek(1:nbr_bdfs,2);
nbr_comb_bdfs           = size(comb_bdfs,1);


% preallocate matrices to store the results
princ_ang_V_potent      = zeros(nbr_comb_bdfs,size(opn_spaces(1).V_potent,2));
princ_ang_V_null        = zeros(nbr_comb_bdfs,size(opn_spaces(1).V_null,2));


% -------------------------------------------------------------------------
% do for all combinations
for i = 1:nbr_comb_bdfs
    
    % compute angles between potent spaces
    princ_ang_V_potent(i,:) = principal_angles( opn_spaces(comb_bdfs(i,1)).V_potent, ...
                            opn_spaces(comb_bdfs(i,2)).V_potent );
    
    % and null spaces
    princ_ang_V_null(i,:) = principal_angles( opn_spaces(comb_bdfs(i,1)).V_null, ...
                            opn_spaces(comb_bdfs(i,2)).V_null );
           
%     % TO DELETE ===========================================================
%     % get angles between manifolds old-style
%     [ang, dim_min_ang]  = comp_hyperplanes_fcn_dim_finding_closest( ...
%                             opn_spaces(comb_bdfs(i,1)).V_potent, ...
%                             opn_spaces(comb_bdfs(i,2)).V_potent, ...
%                             1:size(opn_spaces(comb_bdfs(i,1)).V_potent,2) ); 
%     [ang_rev, dim_min_ang_rev]  = comp_hyperplanes_fcn_dim_finding_closest( ...
%                             opn_spaces(comb_bdfs(i,2)).V_potent, ...
%                             opn_spaces(comb_bdfs(i,1)).V_potent, ...
%                             1:size(opn_spaces(comb_bdfs(i,1)).V_potent,2) ); 
%     % UP TO HERE===========================================================
end


% return variables
comp_opn_spaces.princ_ang_V_potent  = princ_ang_V_potent;
comp_opn_spaces.princ_ang_V_null    = princ_ang_V_null;

% % TO DELETE ===========================================================
% comp_opn_spaces.angles_old_school.angle_12 = ang;
% comp_opn_spaces.angles_old_school.angle_21 = ang_rev;
% comp_opn_spaces.angles_old_school.dim_min_angle_12 = dim_min_ang;
% comp_opn_spaces.angles_old_school.dim_min_angle_21 = dim_min_ang_rev;
% % UP TO HERE===========================================================
