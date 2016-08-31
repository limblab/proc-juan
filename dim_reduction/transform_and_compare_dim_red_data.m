%
% Transform neural data from task 1 ('within') into the PC space of task 2
% ('across'), and compare the projections they yield 
%
%   pc_proj_across_tasks = transform_and_compare_dim_red_data( dim_red_FR, ...
%                               smoothed_FR, labels, neural_chs, within_task, ...
%                               across_task, comp_nbr )
%
% Inputs:
%   dim_red_FR              : cell array with reduced neural data
%   smoothed_FR             : cell array with smoothed FRs
%   labels                  : cell array with a label describing each trial
%   neural_chs              : channels that will be analyzed. If empty, it
%                               will look at all
%   within_task             : task which data will be projected using the
%                               PC matrix for 'across_task' (scalar)
%   across_task             : task which PC decomposition will be used to
%                               project the 'within_task' data
%   comp_nbr                : 1-D vector with the components that will be
%                               compared, or,
%                             2-D vector that established which transformed
%                             dimensions in 'within' will be compared to
%                             which dimensions in 'across'
%   (plot_yn)               : [true] summary plot per component
%
% Note: if one of the eigenvectors has opposite orientation in the two
% tasks, the function already fixes that
%
%

function pc_proj_across_tasks = transform_and_compare_dim_red_data( dim_red_FR, ...
                            smoothed_FR, labels, neural_chs, within_task, ...
                            across_task, comp_nbr, varargin )

                        
% read inputs
plot_yn                 = true;
if nargin == 8
    plot_yn             = varargin{1};
end


% define a couple of parameters
% min_t and max_t for the time domain plots --could be a parameter
t_lims                  = [0 30];
% number of shuffled versions of the projections that will be created to
% assess significance
nbr_corrs_chance        = 100;


% some checks
if length(dim_red_FR) ~= length(smoothed_FR)
    error('length of dim_red_FR and smoothed_FR has to be the same');
end
if length(dim_red_FR) ~= length(labels)
    error('labels has wrong size');
end


% get rid of the first column in smoothed_FR, which is the time
for i = 1:length(smoothed_FR)
    smoothed_FR{i}(:,1)     = [];
end
% get rid of the neural channels we don't want to use (i.e. those that are
% not listed in neural_chs). If neural_chs is empty, it will use all
if ~isempty(neural_chs)
    if length(neural_chs) ~= size(smoothed_FR{1},2)
        for i = 1:length(smoothed_FR)
            smoothed_FR{i}  = smoothed_FR{i}(:,neural_chs+1);
        end
    end
end


% populate a 2-D array that establishes which components in within will be
% comapred to which components in across (1st and 2nd columns)
[nbr_rows, nbr_cols]    = size(comp_nbr);
if min(size(comp_nbr)) > 1
    if nbr_cols ~= 2
        error('comp_nbr can only have 1 or 2 columns');
    else
        comp_nbr_array  = comp_nbr;
    end
else % if comp_nbr is a row vector (=> will comp dim X to dim X)
    if nbr_rows < nbr_cols, comp_nbr = comp_nbr'; end;
    comp_nbr_array      = [comp_nbr, comp_nbr];
end


% -------------------------------------------------------------------------
% some definitions
                        
nbr_comps               = length(comp_nbr);

% get bin width
bin_width_neurons       = mean(diff(dim_red_FR{1}.t));

% % interval for the xcorr, in nbr of bins
% int_xcorr               = 30; 
% % time axis for xcorr
% t_axis_xcorr            = bin_width_neurons*(-int_xcorr:1:int_xcorr);
% % nbr points coherence
% nfft_coh                = 1024/2+1;


% initalize struct for returning results
pc_proj_across_tasks    = struct('within_task',within_task,...
                            'across_task',across_task,...
                            'comp_nbr',comp_nbr_array,...
                            'scores_within',zeros(size(dim_red_FR{within_task}.scores,1),nbr_comps),...
                            'scores_across',zeros(size(dim_red_FR{within_task}.scores,1),nbr_comps),...
                            'R2',zeros(1,nbr_comps),...
                            'R2_shuffled',zeros(nbr_corrs_chance,nbr_comps),...
                            'R2_other_projs',zeros(size(smoothed_FR{1},2)-1,nbr_comps),...
                            'dim_max_other_R2',zeros(1,nbr_comps),...
                            'r',zeros(1,nbr_comps),...
                            'r_shuffled',zeros(nbr_corrs_chance,nbr_comps),...
                            'r_other_projs',zeros(size(smoothed_FR{1},2)-1,nbr_comps),...
                            'dim_max_other_r',zeros(1,nbr_comps));
%                             'vaf',zeros(1,nbr_comps), ...
%                             'vafR',zeros(1,nbr_comps), ...
%                             'CC',zeros(1,nbr_comps), ...
%                             'CCstats',cell(1,nbr_comps));


                      
% -------------------------------------------------------------------------
% do
for i = 1:nbr_comps


    % ---------------------------------------------------------------------
    % Project the data from task 'within task' onto eigenvector 'i' of its
    % manifold
    within_proj             = dim_red_FR{within_task}.scores(:,comp_nbr_array(i,1));
    
    % Project the data from task 'within_task' onto the corresponding
    % eigenvector in the 'across_task' manifold
    across_proj             = smoothed_FR{within_task}*...
                            dim_red_FR{across_task}.w(:,comp_nbr_array(i,2))...
                            - mean(smoothed_FR{within_task})*...
                            dim_red_FR{across_task}.w(:,comp_nbr_array(i,2));

                        
	% ---------------------------------------------------------------------
    % Coherence, cross-correlattion, canonical correlation, and VAF
    
%     % compute cross-correlation
%     xcorr_this_comb         = xcorr( within_proj, across_proj, int_xcorr );
% 
%                             
%     % compute coherence
%     [coh_this_comb, f_coh]  = mscohere( within_proj, across_proj, 20, 16, 1024, 20 );
%     
%     % invert the components that need to be inverted for the plots
%     [~, indx_max_xcorr ]    = max(abs(pc_proj_across_tasks.xcorr(:,i)));
%     if pc_proj_across_tasks.xcorr(indx_max_xcorr,i) < 0
%         inverted_eigenv_this = 1;
%         across_proj         = - across_proj;
%         xcorr_this_comb     = xcorr_this_comb;
%     else
%         inverted_eigenv_this = 0;
%     end
%
%     % compute canonical correlation
%     [~,~,CC_this_comb,~,~,CCstats_this_comb]    = canoncorr( within_proj, across_proj );
%
%     
%     % compute VAF
%     vaf_this_comb           = calc_vaf( within_proj, across_proj ); % limblab fcn
%     vaf_this_combR          = calc_VAF( within_proj, across_proj ); % Raeed's fcn


	% ---------------------------------------------------------------------
    % R2 and correlation between projections
    
    R2_this_comb            = CalculateR2( within_proj, across_proj );
    
    [r_this_comb, P_this_comb] = corrcoef( within_proj, across_proj );
    r_this_comb             = abs(r_this_comb(1,2));
    P_this_comb             = P_this_comb(1,2);
    
    
    % ---------------------------------------------------------------------
    % Control #1: shuffle one of the projections in time, for comparison
    
    
    % * For the correlation
    r_rnd                   = zeros(1,nbr_corrs_chance-1);
    % reverse one of them
    r_rev                   = corrcoef( within_proj, fliplr(across_proj') );
    r_rev                   = abs(r_rev(1,2));
    % and a lot of random perturbations
    for j = 2:nbr_corrs_chance
        r_rnd_this          = corrcoef( within_proj, across_proj(randperm(length(across_proj))));
        r_rnd(j-1)          = abs(r_rnd_this(1,2));
    end
    r_shuffled              = [r_rev, r_rnd];
    
    % * For the R2
    R2_rnd                  = zeros(1,nbr_corrs_chance-1);
    R2_rev                  = CalculateR2( within_proj, fliplr(across_proj')' );
    for j = 2:nbr_corrs_chance
        R2_rnd(j-1)         = CalculateR2( within_proj, across_proj(randperm(length(across_proj))));
    end
    R2_shuffled               = [R2_rev, R2_rnd];
    
    clear r_rnd r_rev R_rnd R_rev
    
    
    % ---------------------------------------------------------------------
    % Control #2: project the neural data onto the other dimensions of the
    % manifold (the 'non corresponding' ones)
    
    indx_ctrl_comp          = setdiff(1:numel(dim_red_FR{1}.chs),comp_nbr_array(i,2));
    
    r_other_projs           = zeros(numel(indx_ctrl_comp),1);
    R2_other_projs          = zeros(numel(indx_ctrl_comp),1);
    for j = 1:numel(indx_ctrl_comp)
        ctrl_across_proj    = smoothed_FR{within_task}*...
                            dim_red_FR{across_task}.w(:,indx_ctrl_comp(j))...
                            - mean(smoothed_FR{within_task})*...
                            dim_red_FR{across_task}.w(:,indx_ctrl_comp(j));
        
        r_this_ctrl_this    = corrcoef( within_proj, ctrl_across_proj );
        r_other_projs(j)    = abs(r_this_ctrl_this(1,2));
        
        R2_other_projs(j)   = CalculateR2( within_proj, ctrl_across_proj );
    end
    
    % find dimension of highest r & R2
    [~, pos_max_r]          = max(r_other_projs);
    dim_max_other_r         = indx_ctrl_comp(pos_max_r);
    [~, pos_max_R2]         = max(R2_other_projs);
    dim_max_other_R2        = indx_ctrl_comp(pos_max_R2);
    
    

    % -----------------------
    % store results in return struct
    pc_proj_across_tasks.scores_within(:,i)     = within_proj;
    pc_proj_across_tasks.scores_across(:,i)     = across_proj;

    pc_proj_across_tasks.R2(i)                  = R2_this_comb;
    pc_proj_across_tasks.R2_shuffled(:,i)       = R2_shuffled;
    pc_proj_across_tasks.R2_other_projs(:,i)    = R2_other_projs;
    pc_proj_across_tasks.dim_max_other_R2(i)    = dim_max_other_R2;
    
    pc_proj_across_tasks.r(i)                   = r_this_comb;
    pc_proj_across_tasks.r_shuffled(:,i)        = r_shuffled;
    pc_proj_across_tasks.r_other_projs(:,i)     = r_other_projs;
    pc_proj_across_tasks.dim_max_other_r(i)     = dim_max_other_r;

%     pc_proj_across_tasks.P_r(i)                 = P_this_comb;
%     pc_proj_across_tasks.vaf(:,i)               = vaf_this_comb;
%     pc_proj_across_tasks.vafR(:,i)              = vaf_this_combR;
%     pc_proj_across_tasks.xcorr(:,i)             = xcorr_this_comb;
%     pc_proj_across_tasks.coh(:,i)               = coh_this_comb;
%     pc_proj_across_tasks.inverted_eigenv(:,i)   = inverted_eigenv_this;
%     if i == 1
%         pc_proj_across_tasks.f_coh              = f_coh;
%     end
%     pc_proj_across_tasks.CC(i)                  = CC_this_comb;
%     pc_proj_across_tasks.CCstats{i}             = CCstats_this_comb;
end



% ------------------------------------------------------------------------
% PLOTS

if plot_yn
    for i = 1:2:nbr_comps
        figure,
        subplot(211),hold on
        plot(dim_red_FR{within_task}.t,[ pc_proj_across_tasks.scores_within(:,i), ...
            pc_proj_across_tasks.scores_across(:,i) ], 'LineWidth',2)
        %plot(dim_red_FR{within_task}.t,dim_red_FR{within_task}.scores(:,comp_nbr(i)),'LineWidth',2)
        legend(labels{within_task},labels{across_task})
        set(gca,'Tickdir','out'),set(gca,'FontSize',14)
        ylabel(['neural comp.' num2str(comp_nbr(i))]),xlim(t_lims)
        subplot(212),hold on
        plot(dim_red_FR{within_task}.t,[ pc_proj_across_tasks.scores_within(:,i+1), ...
            pc_proj_across_tasks.scores_across(:,i+1) ], 'LineWidth',2)
        %plot(dim_red_FR{within_task}.t,dim_red_FR{within_task}.scores(:,comp_nbr(i)),'LineWidth',2)
        legend(labels{within_task},labels{across_task})
        set(gca,'Tickdir','out'),set(gca,'FontSize',14)
        ylabel(['neural comp.' num2str(comp_nbr(i)+1)]),xlabel('time (s)'), xlim(t_lims)

        %         subplot(223)
%         plot(t_axis_xcorr,pc_proj_across_tasks.xcorr(:,i),'LineWidth',2)
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14)
%         ylabel('crosscorrelation'), xlabel('time (s)')
%         subplot(224)
%         plot(f_coh,pc_proj_across_tasks.coh(:,i),'LineWidth',2)
%         set(gca,'Tickdir','out'),set(gca,'FontSize',14)
%         ylabel('coherence'), xlabel('frequency (Hz)'), ylim([0 1])
%         set(gcf,'Colormap',winter)
    end
end