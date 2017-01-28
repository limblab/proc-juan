%
% Plot decoder predictions for data generated with build_neuron_decoder,
% build_dPC_decoder, build_PC_decoder
%
%   function plot_preds_decoder( fd, varargin )
%
% Inputs (opt)      : [default]
%   fd              : decoder data struct
%   (vars)          : ['all'] output variables to plot
%
%

function plot_preds_decoder( fd, varargin )


if nargin == 2
    chs     = varargin{1};
else
    chs     = 1:size(fd.all.y,2);
end

% create time vector
t           = 0:fd.params.bin_size:fd.params.bin_size*(size(fd.all.y,1)-1);


figure,
for v = 1:numel(chs)
    tv      = chs(v);
    subplot(numel(chs),1,v), hold on
    plot(t,fd.all.y_pred(:,tv),'k','linewidth',1.5)
    plot(t,fd.all.y(:,tv),'color',[.5 .5 .5],'linewidth',1.5)
    set(gca,'TickDir','out','FontSize',12), box off
    legend(['R^2=' num2str(fd.all.R2(tv))]), legend boxoff
    ylabel(['Var. ' num2str(tv)])
    if v == numel(chs), xlabel('Time (s)'); end
end