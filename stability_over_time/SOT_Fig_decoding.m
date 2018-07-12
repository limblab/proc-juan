% 
% Plot results of tests of decoding stability. 
%


function SOT_Fig_decoding( dec_results, dec_spike_results )


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

hist_x      = 0:0.05:1.25;

a_to_w      = mn_a./mn_w;
a_to_s_w    = mn_a./mn_s_w;
u_to_s_w    = mn_u./mn_s_w;
s_a_to_s_w  = mn_s_a./mn_s_w;

n_hist      = length(a_to_w);

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


% test statistical comparison
p_norm_a_to_u = ranksum(a_to_s_w,u_to_s_w);
p_norm_a_to_s_a = ranksum(a_to_s_w,s_a_to_s_w);


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
figure, hold on
plot( x, y_w, 'color', cols_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_u, 'color', cols_u,'linewidth', 2)

plot( dd, mn_w,'.','color',cols_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_u,'.','color',cols_u,'markersize',32 )

ylim([0 1]), xlim([0 max(dd)+1]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
legend('within-day','alinged','unaligned','Location','east'), legend boxoff
set(gcf,'color','w'), title('Decoders based on latent activity')


% 2. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY DECODER BASED ON SPIKES
figure, hold on
plot( x, y_w, 'color', cols_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_s_a, 'color', cols_s_a,'linewidth', 2)

plot( dd, mn_w,'.','color',cols_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_s_a,'.','color',cols_s_a,'markersize',32 )

ylim([0 1]), xlim([0 max(dd)+1]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
legend('latent within-day','alinged across','spikes across','Location','east'), legend boxoff
set(gcf,'color','w')



% 3. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY DECODER BASED ON SPIKES
figure, hold on
plot( x, y_s_w, 'color', cols_s_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_s_a, 'color', cols_s_a,'linewidth', 2)

plot( dd, mn_s_w,'.','color',cols_s_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_s_a,'.','color',cols_s_a,'markersize',32 )

ylim([0 1]), xlim([0 max(dd)+1]);
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Days from decoder training'), ylabel('Hand velocity prediction (R^2)')
legend('spikes within-day','alinged across','spikes across','Location','east'), legend boxoff
set(gcf,'color','w')


%% ------------------------------------------------------------------------
% PLOT DISTRIBUTION OF NORMALIZED ALIGNED VS UNALIGNED PREDICTIONS

xh          = hist_x(1:end-1);


% Normalized aligned vs unaligned preds

figure, hold on
h1 = bar(xh,hist_a_to_s_w,'histc');
set(h1,'FaceColor',cols_a)
h2 = bar(xh,hist_u_to_s_w,'histc');
set(h2,'FaceColor',cols_u)
alpha(h2,.5)

tm = max([hist_a_to_s_w,hist_u_to_s_w]);
yln = ceil(tm/5)*5+5;
ys = ceil(tm) + (yln-ceil(tm)) / 2;

plot(mn_a_to_s_w,ys,'.','markersize',32,'color',cols_a)
plot([mn_a_to_s_w-sd_a_to_s_w, mn_a_to_s_w+sd_a_to_s_w],[ys ys],'color',cols_a,'linewidth',3)
plot(mn_u_to_s_w,ys,'.','markersize',32,'color',cols_u)
plot([mn_u_to_s_w-sd_u_to_s_w, mn_u_to_s_w+sd_u_to_s_w],[ys ys],'color',cols_u,'linewidth',3)

text(.1,yln,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,yln-1,['P=' num2str(p_norm_a_to_u,2)],'Fontsize',14)

legend('aligned','unaligned'), legend boxoff
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Normalized prediction accuracy')
ylabel('Session comparisons (%)')
set(gcf,'color','w')


% Normalized aligned vs unaligned preds

figure, hold on
h3 = bar(xh,hist_a_to_s_w,'histc');
set(h3,'FaceColor',cols_a)
h4 = bar(xh,hist_s_a_to_s_w,'histc');
set(h4,'FaceColor',cols_s_a)
alpha(h4,.5)

tm = max([hist_a_to_s_w,hist_s_a_to_s_w]);
yln = ceil(tm/5)*5+5;
ys = ceil(tm) + (yln-ceil(tm)) / 2;

plot(mn_a_to_s_w,ys,'.','markersize',32,'color',cols_a)
plot([mn_a_to_s_w-sd_a_to_s_w, mn_a_to_s_w+sd_a_to_s_w],[ys ys],'color',cols_a,'linewidth',3)
plot(mn_s_a_to_s_w,ys,'.','markersize',32,'color',cols_s_a)
plot([mn_s_a_to_s_w-sd_s_a_to_s_w, mn_s_a_to_s_w+sd_s_a_to_s_w],[ys ys],'color',cols_s_a,'linewidth',3)

text(.1,yln,['n=' num2str(n_hist)],'Fontsize',14)
text(.1,yln-1,['P=' num2str(p_norm_a_to_s_a,2)],'Fontsize',14)

legend('aligned','spikes across'), legend boxoff
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Normalized prediction accuracy')
ylabel('Session comparisons (%)')
set(gcf,'color','w')