% 
% Plot results of tests of decoding stability. 
%


function SOT_Fig_decoding( dec_results, dec_spike_results, params )


%% ------------------------------------------------------------------------
% Get some info

% days between sessions
dd          = dec_results.diff_days;
x           = [0 max(dd)+1];

% stats for the decoders based on latent activity
mn_w        = mean(dec_results.withinR2_m,2)';
mn_a        = mean(dec_results.acrossR2,2)';
mn_u        = mean(dec_results.ctrlR2,2)';


% stats for the decoders based on spikes
mn_s_w      = mean(dec_spike_results.withinR2_m,2)';
mn_s_a      = mean(dec_spike_results.acrossR2,2)';


%% ------------------------------------------------------------------------
% do linear fits of the predictions vs time

lf_w        = polyfit(dd, mn_w, 1);
y_w         = polyval(lf_w,x);

lf_a        = polyfit(dd, mn_a, 1);
y_a         = polyval(lf_a,x);

lf_u        = polyfit(dd, mn_u, 1);
y_u         = polyval(lf_u,x);

lf_s_w      = polyfit(dd, mn_s_w, 1);
y_s_w       = polyval(lf_s_w,x);

lf_s_a      = polyfit(dd, mn_s_a, 1);
y_s_a       = polyval(lf_s_a,x);


%% ------------------------------------------------------------------------
% compute distribution of how well across-day decoders do compared to
% within day decoders

hist_x      = 0:0.05:1.5;

a_to_w      = mn_a./mn_w;
a_to_s_w    = mn_a./mn_s_w;
u_to_s_w    = mn_u./mn_s_w;
s_a_to_s_w  = mn_s_a./mn_s_w;

n_hist      = length(a_to_w);

% histograms
hist_a      = histcounts(mn_a,hist_x)/n_hist*100;
hist_s_w    = histcounts(mn_s_w,hist_x)/n_hist*100; % spikes within
hist_u      = histcounts(mn_u,hist_x)/n_hist*100;

% normalized histograms
hist_a_to_w = histcounts(a_to_w,hist_x)/n_hist*100;
hist_a_to_s_w = histcounts(a_to_s_w,hist_x)/n_hist*100;
hist_u_to_s_w = histcounts(u_to_s_w,hist_x)/n_hist*100;
hist_s_a_to_s_w = histcounts(s_a_to_s_w,hist_x)/n_hist*100;

% stats for the histograms
mn_a_to_w   = mean(a_to_w);
sd_a_to_w   = std(a_to_w);

mn_a_to_s_w = mean(a_to_s_w);
sd_a_to_s_w = std(a_to_s_w);

mn_u_to_s_w = mean(u_to_s_w);
sd_u_to_s_w = std(u_to_s_w);

mn_s_a_to_s_w = mean(s_a_to_s_w);
sd_s_a_to_s_w = std(s_a_to_s_w);

% not normalized
mn_a_hist = mean(mn_a);
mn_u_hist = mean(mn_u);
mn_s_w_hist = mean(mn_s_w);

sd_a_hist = std(mn_a);
sd_u_hist = std(mn_u);
sd_s_w_hist = std(mn_s_w);


% test statistical comparison
p_norm_a_to_u = ranksum(a_to_s_w,u_to_s_w);
p_norm_a_to_s_a = ranksum(a_to_s_w,s_a_to_s_w);

p_a_u = ranksum(mn_a,mn_u);
p_a_w = ranksum(mn_a,mn_w);

%% ------------------------------------------------------------------------
% PLOT PERFORMANCE VS TIME

% colors
cols_w      = [.7 .7 .7];
cols_a      = [0 0 0];
cols_u      = [250,128,114]/255;
cols_s_w    = [175,238,238]/255;
cols_s_a    = [255,99,71]/255;


% 1. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY UNALIGNED
f1 = figure; hold on
plot( x, y_w, 'color', cols_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_u, 'color', cols_u,'linewidth', 2)

plot( dd, mn_w,'.','color',cols_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_u,'.','color',cols_u,'markersize',32 )

ylim([0 1]), xlim([x(1) x(2)]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
legend('within-day','alinged','unaligned','Location','east'), legend boxoff
set(gcf,'color','w'), title('Decoders based on latent activity')


% 2. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY DECODER BASED ON SPIKES
f2 = figure; hold on
plot( x, y_w, 'color', cols_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_s_a, 'color', cols_s_a,'linewidth', 2)

plot( dd, mn_w,'.','color',cols_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_s_a,'.','color',cols_s_a,'markersize',32 )

ylim([0 1]), xlim([x(1) x(2)]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
legend('latent within-day','alinged across','spikes across','Location','east'), legend boxoff
set(gcf,'color','w')



% 3. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY DECODER BASED ON SPIKES
f3 = figure; hold on
plot( x, y_s_w, 'color', cols_s_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_s_a, 'color', cols_s_a,'linewidth', 2)

plot( dd, mn_s_w,'.','color',cols_s_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_s_a,'.','color',cols_s_a,'markersize',32 )

ylim([0 1]), xlim([0 max(dd)+1]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
if prctile(mn_s_a,90) > 0.4
    legend('spikes within-day','aligned across','spikes across','Location','Southwest'), 
else
    legend('spikes within-day','aligned across','spikes across','Location','NorthEast'), 
end
legend boxoff
set(gcf,'color','w')


%% ------------------------------------------------------------------------
% PLOT DISTRIBUTION OF NORMALIZED ALIGNED VS UNALIGNED PREDICTIONS

xh          = hist_x(1:end-1);


%  Aligned vs unaligned preds vs spike preds
f4 = figure; hold on
h0 = bar(xh,hist_u,'histc');
set(h0,'facecolor',cols_u)
alpha(h0,0.5);
h1 = bar(xh,hist_s_w,'histc');
set(h1,'facecolor',cols_s_w)
alpha(h1,0.5);
h2 = bar(xh,hist_a,'histc');
set(h2,'facecolor',cols_a)
alpha(h2,0.5);

tm = max([hist_s_w,hist_a,hist_u]);
yln = ceil(tm/5)*5+5;
ys = ceil(tm) + (yln-ceil(tm)) / 2;

plot(mn_u_hist,ys,'.','markersize',32,'color',cols_u)
plot([mn_u_hist-sd_u_hist, mn_u_hist+sd_u_hist],[ys ys],'color',cols_u,'linewidth',2)
plot(mn_s_w_hist,ys,'.','markersize',32,'color',cols_s_w)
plot([mn_s_w_hist-sd_s_w_hist, mn_s_w_hist+sd_s_w_hist],[ys ys],'color',cols_s_w,'linewidth',2)
plot(mn_a_hist,ys,'.','markersize',32,'color',cols_a)
plot([mn_a_hist-sd_a_hist, mn_a_hist+sd_a_hist],[ys ys],'color',cols_a,'linewidth',2)

text(.1,yln-5,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,yln-8,['P_{a,u}=' num2str(p_a_u,2)],'Fontsize',14)
text(.1,yln-11,['P_{a,w}=' num2str(p_a_w,2)],'Fontsize',14)

legend('unaligned','spikes across','aligned'), legend boxoff
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Prediction accuracy (R^2)')
ylabel('Session comparisons (%)')
set(gcf,'color','w')
xlim([0 1])


% Normalized aligned vs unaligned preds

f5 = figure; hold on
h3 = bar(xh,hist_a_to_s_w,'histc');
set(h3,'FaceColor',cols_a)
h4 = bar(xh,hist_s_a_to_s_w,'histc');
set(h4,'FaceColor',cols_s_a)
alpha(h4,.5)

tm = max([hist_a_to_s_w,hist_s_a_to_s_w]);
yln = ceil(tm/5)*5+5;
ys = ceil(tm) + (yln-ceil(tm)) / 2;

plot(mn_a_to_s_w,ys,'.','markersize',32,'color',cols_a)
plot([mn_a_to_s_w-sd_a_to_s_w, mn_a_to_s_w+sd_a_to_s_w],[ys ys],'color',cols_a,'linewidth',2)
plot(mn_s_a_to_s_w,ys,'.','markersize',32,'color',cols_s_a)
plot([mn_s_a_to_s_w-sd_s_a_to_s_w, mn_s_a_to_s_w+sd_s_a_to_s_w],[ys ys],'color',cols_s_a,'linewidth',2)

text(.1,yln,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,yln-5,['P=' num2str(p_norm_a_to_s_a,2)],'Fontsize',14)

legend('aligned','spikes across'), legend boxoff
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Normalized prediction accuracy')
ylabel('Session comparisons (%)')
set(gcf,'color','w')
xlim([0 1.5])




% SAVE FIG ?
if params.decoder_params.save_fig
    
    % LATENT ACTIVITY ONLY-DECODERS
    fn1 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Decoding_over_time_Latent_activity_only_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f1,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn1 '.fig']));
    saveas(f1,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn1 '.png']));
    saveas(f1,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn1 '.pdf']));
    
    % LATENT ACTIVITY ALIGNED AND WITHIN DAY AND SPIKES
    fn2 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Decoding_over_time_Latent_whitin_day_across_day_Spikes_across_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f2,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn2 '.fig']));
    saveas(f2,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn2 '.png']));
    saveas(f2,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn2 '.pdf']));
    
    % LATENT ACTIVITY ALIGNED AND WITHIN DAY AND SPIKES
    fn3 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Decoding_over_time_Latent_across_day_Spikes_within_day_across_day_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f3,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn3 '.fig']));
    saveas(f3,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn3 '.png']));
    saveas(f3,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn3 '.pdf']));
 
    % ALL DISTRIBUTIONS
    fn4 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Decoding_distribution_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f4,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn4 '.fig']));
    saveas(f4,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn4 '.png']));
    saveas(f4,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn4 '.pdf']));
       
    % DISTRIBUTION LATENT ACTIVITY AND SPIKE DECODERS -NORMALIZED 
    fn5 = [params.monkey '_' params.spiking_inputs{1}(1:end-7) '_Decoding_normalized_distribution_' num2str(length(params.mani_dims)) 'D'];
    
    savefig(f5,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn5 '.fig']));
    saveas(f5,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn5 '.png']));
    saveas(f5,fullfile(params.save_dir,params.spiking_inputs{1}(1:end-7),[fn5 '.pdf']));
    
end
