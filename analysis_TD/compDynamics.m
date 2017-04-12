%
% Do CCA to compare time-varying signals across sessions, targets, etc
%

function cca_info = compDynamics( trial_data, signals, idx1, idx2, varargin )

if nargin >= 1
    chs = varargin{1};
else
    aux = cell2mat(cellfun(@size, {master_td.(signals)},'UniformOutput',false)');
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


% and compare the signals that we want
[A, B, r, U, V] = canoncorr( ts1, ts2 );

% Put return var together
cca_info.cc = r;
cca_info.A  = A;
cca_info.B  = B;
cca_info.U  = U;
cca_info.V  = V;