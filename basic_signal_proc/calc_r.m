% 
% Calculate correlation between two sets of time series
%
%   r = calc_r( s1, s2 )
%
% Inputs:
%   s1          : first set of signals as time x variables
%   s2          : secong set of signals as time x variables
%
% Outputs:
%   r           : correlation 


function r = calc_r( s1, s2 )


num_sigs        = size(s1,2);
r               = zeros(1,num_sigs);

for s = 1:num_sigs
    aux_r       = corrcoef(s1(:,s),s2(:,s));
    r(s)        = aux_r(1,2);
end
