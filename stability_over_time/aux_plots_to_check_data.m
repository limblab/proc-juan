%
% Script to plot raw data
%


%% 
% PLOT SINGLE TRIAL KINEMATICS OF REACHES TO ONE TARGET ACROSS SESSIONS


% Variable to plot
var = 'vel'; % 'vel' 'pos' 'acc'

% Target to plot
target_n        = 1;


time            = 0:master_td(1).bin_size:master_td(1).bin_size*(size(master_td(1).pos,1)-1);

idx_plot        = getTDidx(master_td,{'target_direction',meta.targets(target_n)});
cols_plot       = parula(length(idx_plot));

figure,
subplot(121),hold on,subplot(122),hold on
for i = 1:length(idx_plot)
    subplot(121),plot(time,master_td(idx_plot(i)).(var)(:,1),'color',cols_plot(i,:));
    subplot(122),plot(time,master_td(idx_plot(i)).(var)(:,2),'color',cols_plot(i,:));
end
for i = 1:2
    subplot(1,2,i)
    set(gca,'TickDir','out','FontSize',14), box off
    xlabel('Time (s)')
end
subplot(121),ylabel(['X ' var])
subplot(122),ylabel(['Y ' var])
set(gcf, 'color', [1 1 1]);



%% 
% PLOT MEAN FIRING RATE FOR ONE NEURON PER TARGET FOR ONE SESSION


% Variable to plot
spiking_plot    = pars.spiking_inputs; % 'M1_spikes' 'PMd_spikes' 'S1_spikes'
% What units?
units_plot      = 1:10;
% What session?
sess_n          = 11; 


time = 0:master_td(1).bin_size:master_td(1).bin_size*(size(master_td(1).pos,1)-1);


idx_plot        = getTDidx(master_td,{'date',meta.sessions{1}});
cols_plot       = parula(length(unique([master_td(idx_plot).target_direction]))+1);

td_avg          = trialAverage( master_td(idx_plot),{'target_direction'});


for u = 1:length(units_plot)
    figure,hold on,
    for i = 1:length(td_avg)
        plot(time,td_avg(i).(spiking_plot{1})(:,units_plot(u)),'color',cols_plot(i,:));
    end
    xlabel('Time (s)'); ylabel(['Unit ' num2str(units_plot(u))])
    set(gcf, 'color', [1 1 1])
    set(gca,'TickDir','out','FontSize',14), box off
end




%% 
% PLOT EIGENVALUE DISTRIBUTIONS FOR EACH SESSION


cols = parula(n_sessions+1);

% Populate matrix perc. VAF
units_session = arrayfun(@(x) length(x.eigen), pca_info);
perc_var = nan(n_sessions,max(units_session));

for s = 1:n_sessions
    perc_var(s,1:units_session(s)) = 100*cumsum(pca_info(s).eigen)/sum(pca_info(s).eigen);
end

figure, hold on
for s = 1:n_sessions
    plot(perc_var(s,:),'color',cols(s,:),'linewidth',1)
%    lgn{s} = meta.sessions{s};
end
xlabel('Neural mode'); ylabel(['Cumulative variance explained (%)'])
xlim([0 20]), ylim([0 100])
% legend(lgn,'Location','NorthWest'), legend boxoff
set(gcf, 'color', [1 1 1])
set(gca,'TickDir','out','FontSize',14), box off



%% 
% CLEAR THE VARS WE DON'T NEED

clearvars -except meta master_td pars n_*