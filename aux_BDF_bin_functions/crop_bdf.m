%
% Crop BDFs between certain specific times t_i and t_f
%
%   function cropped_bdf = crop_bdf( bdf_array, varargin )
%
%
% Input (optional):
%   bdf_array               : array of BDFs
%   (t_i)                   : initial time for cropping (s)
%   t_f                     : final time for cropping (s). If passing an
%                               array of BDFs and set to 'min', it will
%                               crop all the BDFs to the minimum common
%                               duration
% Output:
%   cropped_bdf             : cropped BDFs
%
%
% Note: the current version ignores raw analog data and does not crop
% good_kin_data
%


function cropped_bdf = crop_bdf( bdf_array, varargin )

nbr_bdfs                = length(bdf_array);

% read input parameters
if nargin == 2
    % in this case the specified time is the desired end of the BDFs
    if isnumeric(varargin{1})
        t_f             = varargin{1};
        % or the minimum common duration
    elseif ischar(varargin{1})
        if strcmp(varargin{1},'min')
            t_f         = min(arrayfun(@(x)x.meta.duration, bdf_array));
            disp(['minimum common duration is: ' num2str(t_f)]);
        end
    end
    t_i                 = 0; % set to zero for cropping
elseif nargin == 3
    % in this case the specified times are the desired start and end of the
    % BDFs, in this order
    t_i                 = varargin{1};
    t_f                 = varargin{2};
else
    error('The function only takes 2 or 3 input arguments');
end


% ------------------
% Some checks
% check that t_i < t_f
if t_i > t_f
    error('t_i cannot be < t_f');
end;
% check that the duration of each BDF is equal or longer that t_f
for i = 1:nbr_bdfs
    if bdf_array(i).meta.duration < t_f
        error(['duration of BDF(' num2str(i) ') < t_f']);
    end
end


% ------------------
% Crop!
for i = 1:nbr_bdfs
    
    % in META, update duration and add t_i and t_f in FileSepTime
    bdf_array(i).meta.duration  = t_f - t_i;
    bdf_array(i).meta.FileSepTime(1) = t_i;
    bdf_array(i).meta.FileSepTime(2) = t_f;
    bdf_array(i).meta.bdf_info  = [bdf_array(i).meta.bdf_info, ...
        '. Cropped with crop_bdf on ' datestr(now,1)];
    
    % UNITS
    for ii = 1:length(bdf_array(i).units)
        % crop the beginning
        indx_i           = find(bdf_array(i).units(ii).ts<t_i);
        bdf_array(i).units(ii).ts(indx_i) = [];
        bdf_array(i).units(ii).waveforms(indx_i,:) = [];
        bdf_array(i).units(ii).ts = bdf_array(i).units(ii).ts - t_i;
        % crop the end
        indx_f           = find(bdf_array(i).units(ii).ts>=(t_f-t_i));
        bdf_array(i).units(ii).ts(indx_f) = [];
        bdf_array(i).units(ii).waveforms(indx_f,:) = [];
    end
    clear indx*;
    
    % RAW
    % -- the current version ignores raw analog data
    % words
    indx_i               = find(bdf_array(i).raw.words(:,1)<t_i);
    if ~isempty(indx_i)
        bdf_array(i).raw.words(indx_i,:)    = [];
        bdf_array(i).raw.words(:,1)         = bdf_array(i).raw.words(:,1) - t_i;
    end
    indx_f               = find(bdf_array(i).raw.words(:,1)>(t_f-t_i),1,'first');
    bdf_array(i).raw.words(indx_f:end,:)    = [];
    clear indx*;
    
    % encoder data
    if isfield(bdf_array(i).raw,'enc')
        indx_i               = find(bdf_array(i).raw.enc(:,1)<t_i);
        if ~isempty(indx_i)
            bdf_array(i).raw.enc(indx_i,:)  = [];
            bdf_array(i).raw.enc(:,1)       = bdf_array(i).raw.enc(:,1) - t_i;
        end
        indx_f               = find(bdf_array(i).raw.enc(:,1)>(t_f-t_i),1,'first');
        bdf_array(i).raw.enc(indx_f:end,:)  = [];
        clear indx*;
    end
    
    % EMG
    if isfield(bdf_array(i),'emg')
        indx_i               = find(bdf_array(i).emg.data(:,1)<t_i);
        if ~isempty(indx_i)
            bdf_array(i).emg.data(indx_i,:)     = [];
            bdf_array(i).emg.data(:,1)          = bdf_array(i).emg.data(:,1) - t_i;
        end
        indx_f               = find(bdf_array(i).emg.data(:,1)>(t_f-t_i),1,'first');
        bdf_array(i).emg.data(indx_f:end,:)     = [];
        clear indx*;
    end
    
    % FORCE
    if isfield(bdf_array(i),'force')
        if bdf_array(i).meta.lab == 1 % lab 1
            indx_i               = find(bdf_array(i).force.data(:,1)<t_i);
            if ~isempty(indx_i)
                bdf_array(i).force.data(indx_i,:)   = [];
                bdf_array(i).force.data(:,1)        = bdf_array(i).force.data(:,1) - t_i;
            end
            indx_f               = find(bdf_array(i).force.data(:,1)>(t_f-t_i),1,'first');
            bdf_array(i).force.data(indx_f:end,:)   = [];
            clear indx*;
        else % everyone else
            indx_i               = find(bdf_array(i).force(:,1)<t_i);
            if ~isempty(indx_i)
                bdf_array(i).force(indx_i,:)   = [];
                bdf_array(i).force(:,1)        = bdf_array(i).force(:,1) - t_i;
            end
            indx_f               = find(bdf_array(i).force(:,1)>(t_f-t_i),1,'first');
            bdf_array(i).force(indx_f:end,:)   = [];
            clear indx*;
        end
    end
    
    % WORDS
    indx_i               = find(bdf_array(i).words(:,1)<t_i);
    if ~isempty(indx_i)
        bdf_array(i).words(indx_i,:)        = [];
        bdf_array(i).words(:,1)             = bdf_array(i).words(:,1) - t_i;
    end
    indx_f               = find(bdf_array(i).words(:,1)>(t_f-t_i),1,'first');
    bdf_array(i).words(indx_f:end,:)        = [];
    clear indx*;
    
    % DATABURSTS
    for ii = 1:length(bdf_array(i).databursts)
        if bdf_array(i).databursts{ii,1} > t_i
            indx_i = ii-1; break;
        end
    end
    if indx_i > 0 % will not do if all databursts happened after t_i
        bdf_array(i).databursts(1:indx_i,:) = [];
        for ii = 1:length(bdf_array(i).databursts)
            bdf_array(i).databursts{ii,1} = cell2mat(bdf_array(i).databursts(ii,1)) - t_i;
        end
    end
    
    for ii = 1:length(bdf_array(i).databursts)
        if bdf_array(i).databursts{ii,1}(1) > (t_f-t_i)
            indx_f = ii; break;
        end
    end
    if exist('indx_f','var') % if no databurst falls after t_f
        bdf_array(i).databursts(indx_f:end,:) = [];
    end
    clear indx*;
    
    % POS
    if isfield(bdf_array(i),'pos')
        indx_i               = find(bdf_array(i).pos(:,1)<t_i);
        if ~isempty(indx_i)
            bdf_array(i).pos(indx_i,:)          = [];
            bdf_array(i).pos(:,1)               = bdf_array(i).pos(:,1) - t_i;
        end
        indx_f               = find(bdf_array(i).pos(:,1)>(t_f-t_i),1,'first');
        bdf_array(i).pos(indx_f:end,:)          = [];
        clear indx*;
    end
    
    % TARGETS
    if bdf_array(i).meta.lab == 1
        indx_i               = find(bdf_array(i).targets.corners(:,1)<t_i);
        if ~isempty(indx_i)
            bdf_array(i).targets.corners(indx_i,:)  = [];
            bdf_array(i).targets.corners(:,1)       = bdf_array(i).targets.corners(:,1) - t_i;
            bdf_array(i).targets.rotation(indx_i,:) = [];
            bdf_array(i).targets.rotation(:,1)      = bdf_array(i).targets.rotation(:,1) - t_i;
        end
        indx_f               = find(bdf_array(i).targets.corners(:,1)>(t_f-t_i),1,'first');
        bdf_array(i).targets.corners(indx_f:end,:)  = [];
        bdf_array(i).targets.rotation(indx_f:end,:) = [];
        clear indx*;
    elseif bdf_array(i).meta.lab == 3
        indx_i               = find(bdf_array(i).targets.centers(:,1)<t_i);
        if ~isempty(indx_i)
            bdf_array(i).targets.centers(indx_i,:)  = [];
            bdf_array(i).targets.centers(:,1)       = bdf_array(i).targets.centers(:,1) - t_i;
        end
        indx_f               = find(bdf_array(i).targets.centers(:,1)>(t_f-t_i),1,'first');
        bdf_array(i).targets.centers(indx_f:end,:)  = [];
        clear indx*;
    end
    
    % GOOD_KINETIC_DATA
    % ToDo
end

% Return variable
cropped_bdf             = bdf_array;