%
% Crop a dim_red_FR struct and a dim_red_emg struct to the same bins as a
% binned_data struct
%
%   function crop_dim_red_data_based_on_binned( cropped_binned_data, dim_red_FR, varargin)
%
%


function [cropped_dim_red_FR, varargout] = crop_dim_red_data_based_on_binned( ...
                                    cropped_binned_data, dim_red_FR, varargin)


% -------------------------------------------------------------------------
% get optional input parameters


if nargin == 3
    dim_red_emg         = varargin{1}; 
end


% -------------------------------------------------------------------------
% check that all the data have been binned to the same bin size


bin_size_cbd            = min(diff(cropped_binned_data.timeframe));
% round because sometimes the bin size is off by a few us...
bin_size_cbd            = round(bin_size_cbd*100)/100;
bin_size_drfr           = mean(diff(dim_red_FR.t));
if exist('dim_red_emg','var')
    bin_size_emg        = mean(diff(dim_red_emg.t_axis));
end

if abs(bin_size_cbd - bin_size_drfr) > 0.001
    error('bin sizes for binned_data and dim_red_FR are different');
end
if exist('dim_red_emg','var')
    if abs(bin_size_emg - bin_size_cbd) > 0.001
        error('bin sizes for binned_data and dim_red_emg are different');
    end
end


% -------------------------------------------------------------------------
% retrieve indexes for cropping


% 1) for dim_red_FR
indxs_crop_fr           = find(ismember(dim_red_FR.t,cropped_binned_data.timeframe));

% check that we've found all the indexes
if length(cropped_binned_data.timeframe) ~= length(indxs_crop_fr)
    error('some datapoints are missing from the dim_red_FR data struct');
end


% 2) for dim_red_EMG
if exist('dim_red_emg','var')
    indxs_crop_emg      = find(ismember(dim_red_emg.t_axis,cropped_binned_data.timeframe));
    
    % check that we've found all the indexes
    if length(cropped_binned_data.timeframe) ~= length(indxs_crop_emg)
        error('some datapoints are missing from the dim_red_emg data struct');
    end
end


% -------------------------------------------------------------------------
% crop the data and define return variables

cropped_dim_red_FR      = dim_red_FR;
cropped_dim_red_FR.t    = cropped_dim_red_FR.t(indxs_crop_fr);
cropped_dim_red_FR.scores = cropped_dim_red_FR.scores(indxs_crop_fr,:);

if exist('dim_red_emg','var')
    cropped_dim_red_emg = dim_red_emg;
    cropped_dim_red_emg.t_axis  = cropped_dim_red_emg.t_axis(indxs_crop_emg,:);
    cropped_dim_red_emg.scores  = cropped_dim_red_emg.scores(indxs_crop_emg,:);
    
    varargout{1}        = cropped_dim_red_emg;
end
