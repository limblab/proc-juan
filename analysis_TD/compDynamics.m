%
% Do CCA to compare time-varying signals across sessions, targets, etc
%
% This is a quick hack with many things to do:
%   1. Use standard TD functions whenever possible
%   2. maybe more things?

function align_info = compDynamics( trial_data, signals, idx1, idx2, varargin )
% align_method is last input (optional) and can be 'cca' or 'procrustes'
align_method = 'cca';
if length(varargin) >= 1
    chs = varargin{1};
    if length(varargin) > 1
        align_method = varargin{2};
    end
else
    aux = cell2mat(cellfun(@size, {trial_data.(signals)},'UniformOutput',false)');
    chs = 1:min(aux(:,2));
end

% Get the trials that we want
td1     = trial_data(idx1);
td2     = trial_data(idx2);

% keep the variables we want
ts1     = cat(1,td1.(signals));
ts2     = cat(1,td2.(signals));

% and the channels we want
ts1     = ts1(:,chs);
ts2     = ts2(:,chs);


switch lower(align_method)
    case 'cca'
        % and compare the signals that we want
        [A, B, r, U, V] = canoncorr( ts1, ts2 );
        
        % Put return var together
        align_info.cc = r;
        align_info.A  = A;
        align_info.B  = B;
        align_info.U  = U;
        align_info.V  = V;
        
    case 'procrustes'
        
        [d,Z,T] = procrustes(ts1,ts2);
        
        r = zeros(1,length(chs));
        for i = 1:length(chs)
            temp = corrcoef(ts1(:,i),Z(:,i));
            r(i) = temp(1,2);
        end
        
        align_info.cc = r;
        align_info.A = [];
        align_info.B = [];
        align_info.U = ts1;
        align_info.V = Z;
        align_info.T = T;
end




