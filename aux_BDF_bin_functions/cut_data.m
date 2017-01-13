%
% Cut continuous data based on bin numbers that define specific windows
% within a trial. Those bin numbers can be obtained using
% find_cutting_times.m
% 
%   c_data = cut_data( data, cut_b )
%
%
% Inputs (opt)      : [default]
%   data            : data to cut ( bins x variables )
%   cut_b           : bin numbers to cut the data ( trials x 2 ). Columns 1
%                       and 2 contain the start and end of each window
%
%

function c_data = cut_data( data, cut_b )


% find the bins we want to keep
bins_keep           = [];

for w = 1:size(cut_b,1)
    bins_keep       = [bins_keep, cut_b(w,1):cut_b(w,2)];
end

% keep the bins we want
c_data              = data(bins_keep,:);