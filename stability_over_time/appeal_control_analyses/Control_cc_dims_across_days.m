%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  true;

num_iters = 100;

monkey = 'mihili'; % 'chewie','chewie2','jaco','mihili'};
array = 'M1';
avg_dims = 1:4;
mani_dims = 1:10;

load([monkey '_data.mat']);

pars.mani_dims = mani_dims;
pars.align_latent_params.mani_dims = mani_dims;

master_td = smoothSignals(master_td,struct('signals','M1_spikes','calc_fr',true));

dates = unique({master_td.date});

targs = unique([master_td.target_direction]);

master_td_all_trials = [];
for i = 1:length(dates)
    [~,td] = getTDidx(master_td,'date',dates{i});
    
    [td,~] = dimReduce(td,struct('signals','M1_spikes'));
    
    master_td_all_trials = [master_td_all_trials, td];
    
end

% compute cc metric for all days
pars.xval_yn = false;
pars.signals = 'M1_pca';
results = align_latent_activity( master_td_all_trials, pars );



%%

figure;
hold all;

plot_colors = brewermap(max(mani_dims),'RdYlBu');
for i = 1:max(mani_dims)
    temp =  mean(results.cc(:,i),2);
    
    plot(results.diff_days,temp,'.','MarkerSize',10,'Color',plot_colors(i,:));
    
end

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 max(results.diff_days)+1],'YLim',[0 1]);

ylabel('Canonical correlation');
xlabel('Days between sessions');
title([monkey ' - ' array]);

set(gca,'CLim',[1 max(mani_dims)]);
colormap(plot_colors);
colorbar;

if save_figs
    saveas(gcf,fullfile(save_dir,'All CC Dims',[ monkey '_' array '_AllCCDims.fig']));
    saveas(gcf,fullfile(save_dir,'All CC Dims',[ monkey '_' array '_AllCCDims.pdf']));
    saveas(gcf,fullfile(save_dir,'All CC Dims',[ monkey '_' array '_AllCCDims.png']));
end


