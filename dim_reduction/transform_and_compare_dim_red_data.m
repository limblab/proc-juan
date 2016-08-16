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


% min_t and max_t for the plot --could be a parameter
t_lims                  = [0 30];


% some checks
if length(dim_red_FR) ~= length(smoothed_FR)
    error('length of dim_red_FR and smoothed_FR has to be the same');
end
if length(dim_red_FR) ~= length(labels)
    error('labels has wrong size');
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
                            'vaf',zeros(1,nbr_comps), ...
                            'vafR',zeros(1,nbr_comps), ...
                            'R2',zeros(1,nbr_comps));%, ...
%                             'CC',zeros(1,nbr_comps), ...
%                             'CCstats',cell(1,nbr_comps));

% pc_proj_across_tasks    = struct('within_task',within_task,'across_task',across_task,...
%                             'comp_nbr',comp_nbr_array,'scores_within',[],'scores_across',[],...
%                             'xcorr',zeros(length(t_axis_xcorr),nbr_comps),...
%                             't_axis_xcorr',t_axis_xcorr', ...
%                             'coh',zeros(nfft_coh,nbr_comps), ...
%                             'f_coh',zeros(nfft_coh,1), ...
%                             'inverted_eigenv',zeros(1,nbr_comps), ...
%                             'vaf',zeros(1,nbr_comps), ...
%                             'vafR',zeros(1,nbr_comps), ...
%                             'R2',zeros(1,nbr_comps));

                        
% -------------------------------------------------------------------------
% do
for i = 1:nbr_comps
   
    within_proj             = dim_red_FR{within_task}.scores(:,comp_nbr_array(i,1));
    
    % Transform the data from task number 'within_task' using the
    % transformation matrix of task number 'across_task,' and remove mean
    % (for comparison, since matlab does remove the mean)
    % -- in smoothed_FR, add 1 to neural_chs because dim 1 is time
    across_proj             = (smoothed_FR{within_task}(:,neural_chs+1))*...
                            dim_red_FR{across_task}.w(:,comp_nbr_array(i,2))...
                            - mean(smoothed_FR{within_task}(:,neural_chs+1))*...
                            dim_red_FR{across_task}.w(:,comp_nbr_array(i,2));

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
    
    % compute VAF
    vaf_this_comb           = calc_vaf( within_proj, across_proj ); % limblab fcn
    vaf_this_combR          = calc_VAF( within_proj, across_proj ); % Raeed's fcn

    % compute R2 by fitting a linear regression model
%     linear_fit              = fitlm( within_proj, across_proj );
%     R2_this_comb            = linear_fit.Rsquared.Ordinary;
    R2_this_comb            = CalculateR2( within_proj, across_proj );
    
%     % compute canonical correlation
%     [~,~,CC_this_comb,~,~,CCstats_this_comb]    = canoncorr( within_proj, across_proj );
    
    % store results in return struct
    pc_proj_across_tasks.scores_within(:,i)     = within_proj;
    pc_proj_across_tasks.scores_across(:,i)     = across_proj;
%     pc_proj_across_tasks.xcorr(:,i)             = xcorr_this_comb;
%     pc_proj_across_tasks.coh(:,i)               = coh_this_comb;
%     pc_proj_across_tasks.inverted_eigenv(:,i)   = inverted_eigenv_this;
%     if i == 1
%         pc_proj_across_tasks.f_coh              = f_coh;
%     end
    pc_proj_across_tasks.vaf(:,i)               = vaf_this_comb;
    pc_proj_across_tasks.vafR(:,i)              = vaf_this_combR;
    pc_proj_across_tasks.R2(:,i)                = R2_this_comb;
%     pc_proj_across_tasks.CC(i)                  = CC_this_comb;
%     pc_proj_across_tasks.CCstats{i}             = CCstats_this_comb;
end


% % Do Canonical correlation between the projecitons onto within and across
% % PCs
% [A,B,CC_this_comb,U,V,CCstats_this_comb]    = canoncorr( pc_proj_across_tasks.scores_within, ...
%                                                 pc_proj_across_tasks.scores_across );


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