% 
% Compute principal angles between plabes A and B
%
%   [angles, varargout] = principal_angles(A,B)
%
% Inputs:
%   A           : n-by-p matrix that defines a p-dimensional plane in an
%               n-dimensional Euclidean space
%   B           : n-by-q matrix that defines a q-dimensional plane in an
%               n-dimensional Euclidean space
%
% Outputs (opt):
%   angles      : principal angles
%   (C)         : singular values of Qa'*Qb
%
%
%
% Note: largely inspired by this post: http://goo.gl/eMclEw


function [angles, varargout] = principal_angles(A,B)

[Qa,~]          = qr(A,0);
[Qb,~]          = qr(B,0);
C               = svd(Qa'*Qb,0);
angles          = acos(C);

if nargout == 2
    varargout{1} = C;
end