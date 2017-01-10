%
% Align and concatenate binned data in X (e.g., neurons/PCs) and y (e.g.,
% emg, force, pos..) 
%
%   function [X, y] = concatenate_to_regress( input_data, output_data, bins_lag )
%
%
% Inputs:
%   input_data          : 3D matrix with input data (bin number x channel x
%                           trial)
%   output_data         : 3D matrix with output data (bin number x channel
%                           x trial)
%   bins_lag            : lag in nbr of bins. > 0 => input_data leads
%
%

function [X, y] = concatenate_to_regress( input_data, output_data, bins_lag )


% cut 'bins_lag' bins to align the data

 % input data leads
if bins_lag >= 0
    input_data(1:bins_lag,:,:)              = [];
    output_data((end-bins_lag+1):end,:,:)   = [];
% output data leads
else
    output_data(1:bins_lag,:,:)             = [];
    input_data((end-bins_lag+1):end,:,:)    = [];
end


% rearrange data matrices to have time vs variable 2D matrices
input_data          = permute(input_data,[1 3 2]);
output_data         = permute(output_data,[1 3 2]);

% and create return vars
X                   = reshape(input_data,[],size(input_data,3));
y                   = reshape(output_data,[],size(output_data,3));
