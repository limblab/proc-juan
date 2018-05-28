% 
% Create a histogram with the % variance accounted for by the manifold
%


if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


mani_dim            = 12;


% -------------------------------------------------------------------------
% Retrieve data into task x dimension matrix 'var_expl'

var_expl            = [];

for d = 1:length(datasets)
    eigenv          = cell2mat( cellfun( @(x) x.eigen, datasets{d}.dim_red_FR, 'UniformOutput', false ) );
    for t = 1:size(eigenv,2)
        cum_eigenv(:,t) = cumsum(eigenv(1:mani_dim,t))/sum(eigenv(:,t));
    end
    var_expl        = [var_expl; cum_eigenv'];
    clear cum_eigenv;
end

% make percentage
var_expl            = var_expl*100;


% -------------------------------------------------------------------------
% Compute histogram
x_hist              = 0:5:100;
c_hist              = histcounts(var_expl(:,end),x_hist);

% Turn counts to percentage
c_hist              = c_hist/sum(c_hist)*100;

% -------------------------------------------------------------------------
% Compute stats
m_vaf               = mean(var_expl(:,end));
sd_vaf              = std(var_expl(:,end));

y_stats             = max(c_hist) + 2.5;


figure,hold on
h_ve = bar(x_hist(1:end-1),c_hist,'histc');
plot([m_vaf-sd_vaf, m_vaf+sd_vaf],[y_stats y_stats],'linewidth',1.5,'color',[.7 .7 .7])
plot(m_vaf,y_stats ,'.','color',[.7 .7 .7],'markersize',20)
set(h_ve,'FaceColor',[.7 .7 .7],'EdgeColor',[.3 .3 .3])
set(gca,'FontSize',14,'TickDir','out'), box off
xlabel('Neural variance expl. (%)'),ylabel('Tasks (%)')
text(10,max(c_hist)-1,['n = ' num2str(size(var_expl,1))],'FontSize',14)



% Clear vars
clearvars -except *_results *_params datasets manifold_dim