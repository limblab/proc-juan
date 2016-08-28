% 
% Load a series of BDFs, process them to run "manifold" analyses, and store
% them into a file with all the necessary info. The data will be stored in
% the same path were the data are, in a file called
% "all_manifold_datasets.mat" 
%
% The preprocessing implies a sequence of: 1) binning the data; 2)
% smoothing the firing rates; 3) doing PCA/FA on the processed spike trains
%
%   function batch_preprocess_dim_red_data( varargin )
%
%
% Inputs (opt)      : [default]
%   (path)          : ['.'] path to the files. Files need to be called
%                       'all_tasks_monkey_date.mat'
%   (params)        : struct with preprocessing params, of type
%                       'batch_preprocess_dim_red_data_params_defaults'
%s
%
% function batch_preprocess_dim_red_data( path )
% function batch_preprocess_dim_red_data( path, params )
%
%


function datasets = batch_preprocess_dim_red_data( varargin )


% -------------------------------------------------------------------------
% read input params

if nargin >= 1
    path            = varargin{1};
else
    path            = pwd;
end

if nargin == 2
    params          = batch_preprocess_dim_red_data_params_defaults(...
                        varargin{2});
else
    params          = batch_preprocess_dim_red_data_params_defaults();
end


% -------------------------------------------------------------------------
% find all "all_tasks_...." files in the folder

files_in_path       = dir(path);
files_indx          = arrayfun( @(x) strncmp(x.name,'all_tasks',9), ...
                        files_in_path, 'UniformOutput', false );
files_indx          = cell2mat(files_indx);
files_indx          = find(files_indx);


% -------------------------------------------------------------------------
% load them, and bin the files, smooth the firing rates and do
% dimensionality reduction. Store meta data as well

datasets            = cell(numel(files_indx),1);

for i = 1:numel(files_indx)
   
    % load the file
    file_name       = files_in_path(files_indx(i)).name;
    load(file_name);
    
    disp(' ');
    disp(['preprocessing file ' file_name]);
    
    % store chosen EMGs and neural channels
    datasets{i}.chosen_emgs     = chosen_emgs;
    datasets{i}.neural_chs      = neural_chs;
    % store task labels, monkey and data
    datasets{i}.labels          = labels;
    underscores_indx            = find(file_name=='_');
    datasets{i}.monkey          = file_name(11:underscores_indx(end)-1);
    datasets{i}.date            = file_name(underscores_indx(end)+1:end);
    clear underscores_indx
    
    % bin the data and smooth the FRs
    %[datasets{i}.smoothed_FR, datasets{i}.binned_data] = gaussian_smoothing( ...
    for ii = 1:numel(cbdf)
        [sfr, bd]               = gaussian_smoothing( ...
                            cbdf(ii), params.transform, ...
                            params.bin_size, params.kernel_SD );
        datasets{i}.smoothed_FR{ii} = sfr;
        datasets{i}.binned_data{ii} = bd;
    end
    
    % do dimensionality reduction on the smoothed FRs
    discard_neurons             = setdiff( 1:numel(cbdf(1).units), neural_chs );
    for ii = 1:numel(cbdf)
        % decide what data to do PCA on --the continous data or the
        % trial-related part of the data only
        switch params.data_pca
            case 'all'
                datasets{i}.dim_red_FR{ii} = dim_reduction( datasets{i}.smoothed_FR{ii}, ...
                            params.dim_red_method, discard_neurons );
            case 'trial-related'
                % split the data in single trials
                % -- note that get-single_trial_data will already exclude
                % the unwanted neural and EMG
                stdata          = get_single_trial_data( datasets{i}.binned_data{ii}, ...
                            labels(ii), neural_chs, chosen_emgs, params.norm_trial_data, ...
                            params.w_i, params.w_f );
                % concatenate the firing rates to do PCA
                stdata          = concatenate_single_trials( stdata );
                % get concatenated smoothed FRs, and add time as first
                % variable --for compatibility with the PCA code
                sfr             = [stdata.target{end}.conc_t', ...
                            stdata.target{end}.neural_data.conc_smoothed_fr];
                % Do PCA on the neural data
                datasets{i}.dim_red_FR{ii} = dim_reduction( sfr, ...
                            params.dim_red_method, [], false, false );
        end
        % add what data were used
        datasets{i}.dim_red_FR{ii}.data = params.data_pca;
    end
    
    % clean up
    clear cbdf bd sfr chosen_emgs neural_chs discard_neurons labels
    
    disp(' ');
end


% save all the datasets
save('all_manifold_datasets.mat','datasets','-v7.3');
disp(['data saved in ' pwd filesep 'all_manifold_datasets.mat']);