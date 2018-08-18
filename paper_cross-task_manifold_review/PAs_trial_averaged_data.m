%
% Signif PAs trial averaged
%

% Load data if it isn't in the workspace
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end


% -------------------------------------------------------------------------
% manifold dimensionality
mani_dims = 12;

% Preallocate stuff
t_session = [];
PAs_single_trials = [];
PAs_trial_avg = [];

ptr = 1;


% DO
for d = 1:length(datasets)
    
    comb_tasks = nchoosek( 1:length(datasets{d}.labels), 2);
    n_comb_tasks = size(comb_tasks,1);
    
    
    % do for all task combinations
    for c = 1:n_comb_tasks
       
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
        
        % -----------------------------------------------------------------
        % For the single trial data
%         W1st = datasets{d}.dim_red_FR{t1}.w(:,1:mani_dims);
%         W2st = datasets{d}.dim_red_FR{t2}.w(:,1:mani_dims);
        
        FR1st = datasets{d}.stdata{t1}.target{end}.neural_data.conc_smoothed_fr;
        FR2st = datasets{d}.stdata{t2}.target{end}.neural_data.conc_smoothed_fr;

        W1st = pca(FR1st);
        W2st = pca(FR2st);
        
        W1st = W1st(:,1:mani_dims);
        W2st = W2st(:,1:mani_dims);
        
        PAst = principal_angles(W1st,W2st);
        
        % -----------------------------------------------------------------
        % for the trial averaged data
         
        FR1ta = datasets{d}.stdata{t1}.target{end}.neural_data.smoothed_fr_mn;
        FR2ta = datasets{d}.stdata{t2}.target{end}.neural_data.smoothed_fr_mn;
        
        W1ta = pca(FR1ta);
        W2ta = pca(FR2ta);
        
        W1ta = W1ta(:,1:mani_dims);
        W2ta = W2ta(:,1:mani_dims);
        
        PAta = principal_angles(W1ta,W2ta);
        
        % -----------------------------------------------------------------
        % save results
        PAs_single_trials = [PAs_single_trials; PAst'];
        PAs_trial_avg = [PAs_trial_avg; PAta'];
        t_session = [t_session, d];
    end
end



% -------------------------------------------------------------------------
% Summary calculations


ex_wrist = 8;
ex_grasp = 11;

% create summary vectors
conc_PAs_single_trial = rad2deg(reshape(PAs_single_trials,1,[]));
conc_PAs_trial_avg = rad2deg(reshape(PAs_trial_avg,1,[]));


% Linear fit scatter plot
lfit = polyfit(conc_PAs_single_trial,conc_PAs_trial_avg,1);
xfit = [0 90];
yfit = polyval(lfit,xfit);


% data distributions
bin_hist = 5;
x_hist = 0:bin_hist:(90+bin_hist);
y_hist_st = histcounts(conc_PAs_single_trial,x_hist)/numel(conc_PAs_single_trial)*100;
y_hist_ta = histcounts(conc_PAs_trial_avg,x_hist)/numel(conc_PAs_trial_avg)*100;


% are the distributions statistically different?
p_dist = ranksum( conc_PAs_single_trial, conc_PAs_trial_avg );


% Scatter plot
figure, hold on
plot(xfit,yfit,'k','linewidth',1.5)
plot([0 90],[0 90],'color',[.6 .6 .6],'linewidth',1.5)
plot(conc_PAs_single_trial,conc_PAs_trial_avg,'.k','markersize',20)
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlabel('Single trial-based principal angles (deg)')
ylabel('PSTH-based principal angles (deg)')
legend('data fit','unity line','location','northwest'), legend boxoff

% Histograms
figure, hold on
b1 = bar(x_hist(1:end-1),y_hist_st,'histc');
set(b1,'FaceColor','b'); alpha(b1,0.6);
b2 = bar(x_hist(1:end-1),y_hist_ta,'histc');
set(b2,'FaceColor',[.6 .6 .6]); alpha(b2,0.6);
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlim([0 90])
ylabel('Principal angle comparisons (%)')
xlabel('Principal angle (deg)')

% Example sessions
t_PAs = rad2deg( PAs_single_trials(t_session == ex_wrist,:));
t_PAs_ta = rad2deg( PAs_trial_avg(t_session == ex_wrist,:));
cols = parula(size(t_PAs,1));
% Example wrist
figure, subplot(122),hold on
for ep = 1:size(t_PAs,1)
    plot( t_PAs(ep,:),'linewidth', 1, 'color', cols(ep,:) )
    plot( t_PAs_ta(ep,:),'-.','linewidth', 1, 'color', cols(ep,:) )
end
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',14,'TickDir','out'); box off
xlabel('Neural mode')
% Example grasp
subplot(121), hold on
t_PAs = rad2deg( PAs_single_trials(t_session == ex_grasp,:));
t_PAs_ta = rad2deg( PAs_trial_avg(t_session == ex_grasp,:));
cols = [0 0 0];
plot( t_PAs, 'linewidth',1,'color','k' )
plot( t_PAs_ta, '-.','linewidth',1,'color','k' )
set(gca,'FontSize',14,'TickDir','out'); box off
xlabel('Neural mode')
ylabel('Principal angle (deg)')
