% 
% Test decoding stability
%


% retrieve some time info
dd          = dec_results.diff_days;
x           = [0 max(dd)+1];

% decoders based on latent activity
mn_w        = mean(dec_results.withinR2_m,2)';
mn_a        = mean(dec_results.acrossR2,2)';
mn_c        = mean(dec_results.ctrlR2,2)';

% decoders based on spikes
mn_s_w      = mean(dec_spike_results.withinR2_m,2)';
mn_s_a      = mean(dec_spike_results.acrossR2,2)';


% linear fits
lf_w        = polyfit(dd, mn_w, 1);
y_w         = polyval(lf_w,x);

lf_a        = polyfit(dd, mn_a, 1);
y_a         = polyval(lf_a,x);

lf_c        = polyfit(dd, mn_c, 1);
y_c         = polyval(lf_c,x);

lf_s_w      = polyfit(dd, mn_s_w, 1);
y_s_w       = polyval(lf_s_w,x);

lf_s_a      = polyfit(dd, mn_s_a, 1);
y_s_a       = polyval(lf_s_a,x);


% colors
cols_w      = [.7 .7 .7];
cols_a      = [0 0 0];
cols_c      = [250,128,114]/255;
cols_s_w    = [175,238,238]/255;
cols_s_a    = [255,99,71]/255;


% 1. DECODERS BASED ON LATENT ACTIVITY: WITHIN-DAY, ACROSS-DAY ALIGNED, AND
% ACROSS-DAY UNALIGNED
figure, hold on
plot( x, y_w, 'color', cols_w,'linewidth', 2)
plot( x, y_a, 'color', cols_a,'linewidth', 2)
plot( x, y_c, 'color', cols_c,'linewidth', 2)

plot( dd, mn_w,'.','color',cols_w,'markersize',32 )
plot( dd, mn_a,'.','color',cols_a,'markersize',32 )
plot( dd, mn_c,'.','color',cols_c,'markersize',32 )

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