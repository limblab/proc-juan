%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  true;

num_iters = 100;

monkey = 'chewie'; % 'chewie','chewie2','jaco','mihili'};
array = 'M1';
avg_dims = 1:4;
mani_dims = 1:10;


load([monkey '_data.mat']);


pars.mani_dims = mani_dims;
pars.align_latent_params.mani_dims = mani_dims;


master_td = smoothSignals(master_td,struct('signals','M1_spikes','calc_fr',true));

dates = unique({master_td.date});

targs = unique([master_td.target_direction]);

[cc,r] = deal(cell(1,num_iters));
for iter = 1:num_iters
    disp(iter);
    
    master_td_all_trials = [];
    for i = 1:length(dates)
        [~,td] = getTDidx(master_td,'date',dates{i});
        
        [td,~] = dimReduce(td,struct('signals','M1_spikes'));
        
        for t = 1:length(targs)
            targ_idx = getTDidx(td,'target_direction',targs(t));
            idx = randperm(length(targ_idx));
            td(targ_idx) = td(targ_idx(idx));
        end
        
        master_td_all_trials = [master_td_all_trials, td];
        
    end
    
    % compute cc metric for all days
    pars.xval_yn = false;
    pars.signals = 'M1_pca';
    results = align_latent_activity( master_td_all_trials, pars );
    
    
    cc{iter} = results.cc;
    r{iter} = results.r;
    
end

diff_days = results.diff_days;

% save it
save(fullfile(save_dir,'Shuffle Trials',[ monkey '_shuffle_trial_order_results.mat'],'cc','r','diff_days'));

%%


across = zeros(length(diff_days),length(cc));
for i = 1:length(cc)
    across(:,i) =  mean(cc{i}(:,avg_dims),2);
end


figure;
hold all;

for iter = 1:size(across,2)
    m = mean(across(:,iter));
    s = std(across(:,iter));
    plot(iter,m,'ko','LineWidth',1);
    plot(iter*[1 1],[m-s, m+s],'k-','LineWidth',2);
end

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 size(across,2)+1],'YLim',[0 1]);

ylabel('Canonical correlation');
xlabel('Iteration Number');
title([monkey ' - ' array]);


if save_figs
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials.fig']));
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials.pdf']));
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials.png']));
end

%%
figure;
hold all;
[n,x] = hist(reshape(across,numel(across),1),0:0.01:1);
h = bar(x,n./sum(n)*100,'hist');
set(h,'FaceColor','k','FaceAlpha',0.3);

set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 1]);

xlabel('Canonical correlation');
ylabel('Density (%)');
title([monkey ' - ' array]);


if save_figs
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials_Hist.fig']));
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials_Hist.pdf']));
    saveas(gcf,fullfile(save_dir,'Shuffle Trial Order',[ monkey '_' array '_ShuffleTrials_Hist.png']));
end



