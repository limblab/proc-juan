
function corr_info = corrDynamics( trial_data, signals, idx1, idx2, varargin )

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

r       = zeros(1,length(chs));

% and compare the signals that we want
for c = 1:length(chs)
    auxcorr = corrcoef( ts1(:,c), ts2(:,c) );
    r(c) = auxcorr(1,2);
end

% Put return var together
corr_info.r = r;