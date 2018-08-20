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


% check that the parallel pool is running, otherwise start it
gcp;


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
        datasets{i}.binned_data{ii} = bd;
    end
    
    % do dimensionality reduction on the smoothed FRs
    if min(size(neural_chs)) == 1 % for threshold crossings
        discard_neurons         = setdiff( 1:numel(cbdf(1).units), neural_chs );
    else % for sorted units
        these_neurons           = cell2mat(arrayfun(@(x) x.id, cbdf(1).units, 'uniformoutput',false)');
        discard_neurons         = [];
        for n = 1:max(size(these_neurons)) 
            % check that a unit in the desired channel and with the desired
            % id is present in the BDF
            if max( sum( these_neurons(n,:) == neural_chs, 2) ) < 2
            
                discard_neurons     = [discard_neurons; n];
            end
        end
    end
    for ii = 1:numel(cbdf)
        % decide what data to do PCA on --the continous data or the
        % trial-related part of the data only
        switch params.data_pca
            
            case 'all'
                datasets{i}.dim_red_FR{ii} = dim_reduction( datasets{i}.smoothed_FR{ii}, ...
                            params.dim_red_method, discard_neurons );
                if params.dim_red_emg ~= 'none'
                    warning('Dim reduction of the EMG not implemented for continuous data');
                end
                
            case 'trial-related'
                
                % ----------------------------------
                % 1) split the data in single trials
                % -- note that get-single_trial_data will already exclude
                % the unwanted neural and EMG channels
                stdata          = get_single_trial_data( datasets{i}.binned_data{ii}, ...
                            labels(ii), neural_chs, chosen_emgs, params.norm_trial_data, ...
                            params.w_i, params.w_f );
                        
                % ----------------------------------
                % 2) PCA of the neural data
                
                % concatenate the firing rates to do PCA
                stdata          = concatenate_single_trials( stdata );
                % get concatenated smoothed FRs, and add time as first
                % variable --for compatibility with the PCA code
                sfr             = [stdata.target{end}.conc_t', ...
                            stdata.target{end}.neural_data.conc_smoothed_fr];
                % Do PCA on the neural data
                datasets{i}.dim_red_FR{ii} = dim_reduction( sfr, ...
                            params.dim_red_method, [], false, false );
                       
                % ----------------------------------
                % 3) Dimensionality reduction of the EMG data
                
                switch params.dim_red_emg
                    case 'none'
                        % do nothing
                    otherwise
                        % do dim reduction if haven't specified to do it
                        % across tasks
                        if ~params.dim_red_emg_across_taks
                            % get concatenated EMG data
                            emg =[stdata.target{end}.conc_t', ...
                                    stdata.target{end}.emg_data.conc_emg];
                            % do PCA/NMF of the EMGs
                            datasets{i}.dim_red_emg{ii} = dim_reduction_muscles( ...
                                emg, params.dim_red_emg, [], [], params.emg_factors );
                        end
                end
                % ----------------------------------
                % 4) Rearrange results
                
                % add dimensionally reduced data to single trial data
                % struct
                if ~isfield(datasets{i},'dim_red_emg')
                    % add the dim_red_FR (neurons) only
                    stdata      = add_dim_red_data_to_single_trial_data( ...
                            stdata, datasets{i}.dim_red_FR{ii} ); 
                elseif ~params.dim_red_emg_across_taks
                    % add the dim_red_FR & dim_red_emg
                    stdata      = add_dim_red_data_to_single_trial_data( ...
                            stdata, datasets{i}.dim_red_FR{ii}, ...
                            datasets{i}.dim_red_emg{ii} ); 
                else
                    % only add dim_red_FR, if want to do PCA/NMF of the
                    % EMGs for all the tasks
                    stdata      = add_dim_red_data_to_single_trial_data( ...
                            stdata, datasets{i}.dim_red_FR{ii} ); 
                end
                
                % add dim_red data to single trial data struct
                % add single trial data to the dataset struct
                datasets{i}.stdata{ii} = stdata;
                
        end
        % add what data were used
        datasets{i}.dim_red_FR{ii}.data = params.data_pca;
        
        % add smoothed_FRs -- it will be the cropped one, if the data were
        % cropped
        datasets{i}.smoothed_FR{ii} = sfr;
    end
    
    
    % do dim reduction of the EMGs across tasks, if specified
    if ~strcmp(params.dim_red_emg,'none') && params.dim_red_emg_across_taks
        tmp_dtst    = add_dim_red_emg_pooled_across_tasks( datasets{i}, params );
        datasets{i}   = tmp_dtst;
    end
    % dim_red_emg = dim_reduction_muscles_pooling_tasks( data_set, 'nmf', [], 4, 0 )
    
    
    % clean up
    clear cbdf bd sfr chosen_emgs neural_chs discard_neurons labels
    
    disp(' ');
end


% save all the datasets
if params.save_data
    save('all_manifold_datasets.mat','datasets','-v7.3');
    disp(['data saved in ' pwd filesep 'all_manifold_datasets.mat']);
end