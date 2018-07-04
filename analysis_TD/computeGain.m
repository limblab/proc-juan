%
% Get gain between two time-varying signals by computing the slope of a
% linear fit (sig2/sig1). If the input signals are M by N matrices, it will
% compute N gains
%

function gain = computeGain( sig1, sig2 )


if size(sig1) ~= size(sig2), error('computeGain: Input signals need to have the same size'); end

gain = zeros(1,size(sig1,2));

for s = 1:size(sig1,2)
    lfit = polyfit(sig1(:,s),sig2(:,s),1);
    gain(s) = lfit(1);
end
    
end