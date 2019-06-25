%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  false;

which_array = 'M1_spikes';

avg_dims = 1:4;
noise_amount = 0.1;

monkeys = {'chewie','chewie2','jaco','mihili'};

[within_vals_rr, ...
    within_vals_tt, ...
    within_align_rr, ...
    within_align_tt, ...
    across_unalign_tr, ...
    across_unalign_tt, ...
    across_align_tr, ...
    across_align_tt, ...
    across_align_rr, ...
    model_align, ...
    model_cid_align, ...
    all_r] = deal([]);
for iMonk = 1:length(monkeys)
    monkey = monkeys{iMonk};
    
    load([monkey '_data.mat']);
    
    pars.mani_dims = 1:10;
    pars.align_latent_params.mani_dims = 1:10;
    
    master_td = smoothSignals(master_td,struct('signals',which_array,'width',pars.kernel_SD,'calc_rate',true));
    master_td_all_trials = smoothSignals(master_td_all_trials,struct('signals',which_array,'width',pars.kernel_SD,'calc_rate',true));
    
    
    % fit tuning curves for all cells to velocity
    %   TCs are fit on whole reach, so they have the "movement information" but
    %   aren't tied to the dynamics
    dates = unique({master_td.date});
    
    td_out = cell(length(dates),2);
    master_td_new = [];
    master_td_new_tune = [];
    for iDate = 1:length(dates)
        [~,td] = getTDidx(master_td,'date',dates{iDate});
        
        
        % td = dupeAndShift(td,'vel',-3);
        
        td_temp = trimTD(td,'start',{'end',-floor(size(td(1).(which_array),1)/2)});
        fr = getSig(td_temp,which_array);
        v = getSig(td_temp,'vel');
        s = sqrt(v(:,1).^2 + v(:,2).^2);
        
        th = atan2(v(:,2),v(:,1));
        X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
        
        tc1 = zeros(size(fr,2),4);
        r = zeros(size(fr,2),1);
        for unit = 1:size(fr,2)
            [b,~,~,~,temp] = regress(fr(:,unit),X);
            r(unit) = temp(1);
            tc1(unit,:) = b;
        end
        
        
        td_temp = trimTD(td,{'start',floor(size(td(1).(which_array),1)/2)+1},'end');
        fr = getSig(td_temp,which_array);
        v = getSig(td_temp,'vel');
        s = sqrt(v(:,1).^2 + v(:,2).^2);
        
        th = atan2(v(:,2),v(:,1));
        X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
        
        tc2 = zeros(size(fr,2),4);
        r = zeros(size(fr,2),1);
        for unit = 1:size(fr,2)
            [b,~,~,~,temp] = regress(fr(:,unit),X);
            r(unit) = temp(1);
            tc2(unit,:) = b;
        end
        
        fr = getSig(td,which_array);
        v = getSig(td,'vel');
        s = sqrt(v(:,1).^2 + v(:,2).^2);
        
        th = atan2(v(:,2),v(:,1));
        % project kinematics through TC to get fake firing rates for all reaches
        [fr_model,fr_model_cid] = deal(zeros(size(fr)));
        tc_switch = repmat(cat(1,zeros(floor(size(td(1).(which_array),1)/2),1),ones(floor(size(td(1).(which_array),1)/2),1)),length(td),1);
        for unit = 1:size(fr,2)
            fr_model(:,unit) = ~tc_switch.*(tc1(unit,1) + ...
                tc1(unit,4)*s + ...
                tc1(unit,2).*s.*sin(th) + ...
                tc1(unit,3).*s.*cos(th) + ...
                -noise_amount*abs(tc1(unit,1))+normrnd(0,2*noise_amount*abs(tc1(unit,1)),size(fr,1),1)) + ...
                tc_switch.*(tc2(unit,1) + ...
                tc2(unit,4)*s + ...
                tc2(unit,2).*s.*sin(th) + ...
                tc2(unit,3).*s.*cos(th) + ...
                -noise_amount*abs(tc2(unit,1))+normrnd(0,2*noise_amount*abs(tc2(unit,1)),size(fr,1),1)) ;
        end
        

        % add modeled FR to TD
        td2 = td;
        count = 0;
        for trial = 1:length(td)
            td2(trial).(which_array) = fr_model(count+1:count+size(td2(trial).(which_array),1),:);
            count = count + size(td2(trial).(which_array),1);
        end
        
        td = dimReduce(td,which_array);
        td2 = dimReduce(td2,which_array);
        
        if any(any(isnan(getSig(td2,'M1_pca'))))
            keyboard
        end
        
        all_r = [all_r; r];
        
        % subsample for within-day
        idx = getTDidx(master_td_all_trials,'date',dates{iDate});
        
        master_td_new = [master_td_new, td];
        master_td_new_tune = [master_td_new_tune, td2];
        
        td_out{iDate,1} = td;
        td_out{iDate,2} = td2;
    end
    
    
    % do within day
    within_align_rr = [within_align_rr, align_latent_activity_within_day(master_td_new, pars.align_latent_params )];
    within_align_tt = [within_align_tt, align_latent_activity_within_day(master_td_new_tune, pars.align_latent_params )];
    
    
    % now alignment across days for tuning model
    for iDate1 = 1:size(td_out,1)
        td = td_out{iDate1,2};
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                td2 = td_out{iDate2,2};
                across_align_tt = [across_align_tt, compDynamics( [td,td2], [which_array(1:end-7) '_pca'], 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
            end
        end
    end
    
    % now alignment across days between tuning modle and real
    for iDate1 = 1:size(td_out,1)
        td = td_out{iDate1,1};
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                td2 = td_out{iDate2,2};
                across_align_tr = [across_align_tr, compDynamics( [td,td2], [which_array(1:end-7) '_pca'], 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
            end
        end
    end
    
     % now alignment across days between real and real
    for iDate1 = 1:size(td_out,1)
        td = td_out{iDate1,1};
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                td2 = td_out{iDate2,1};
                across_align_rr = [across_align_rr, compDynamics( [td,td2], [which_array(1:end-7) '_pca'], 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
            end
        end
    end
    
    % now unaligned correlation across days between tuning model and real
    for iDate1 = 1:size(td_out,1)
        td = td_out{iDate1,1};
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                td2 = td_out{iDate2,2};
                across_unalign_tr = [across_unalign_tr, corrDynamics( [td,td2], [which_array(1:end-7) '_pca'], 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
            end
        end
    end
    
    % now unaligned correlation across days between tuning model
    for iDate1 = 1:size(td_out,1)
        td = td_out{iDate1,2};
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                td2 = td_out{iDate2,2};
                across_unalign_tt = [across_unalign_tt, corrDynamics( [td,td2], [which_array(1:end-7) '_pca'], 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
            end
        end
    end
    
    % pull out the values for normalization
    for iDate1 = 1:size(td_out,1)
        for iDate2 = 1:size(td_out,1)
            if iDate1 ~= iDate2
                within_vals_rr = [within_vals_rr; within_align_rr(iDate2).cc_m];
                within_vals_tt = [within_vals_tt; within_align_tt(iDate2).cc_m];
            end
        end
    end
    
end

%%
figure('Position',[100 100 1000 400]);
% tuning-to-tuning within
[vals_within_tt] = deal(zeros(1,length(within_align_rr)));
for iDate = 1:length(within_align_tt)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align_tt(iDate).aligned_info.cc});
   vals_within_tt(iDate) = mean(temp);
end
% real data within
[val_within_rr] = deal(zeros(1,length(within_align_rr)));
for iDate = 1:length(within_align_rr)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align_rr(iDate).aligned_info.cc});
    val_within_rr(iDate) = mean(temp);
end
% tune to tune across
[val_across_tt] = deal(zeros(1,length(across_align_tt)));
for iDate = 1:length(across_align_tt)
    temp = cellfun(@(x) mean(x(avg_dims)),{across_align_tt(iDate).cc});
    val_across_tt(iDate) = mean(temp);
end
% tune to real across
[val_across_tr] = deal(zeros(1,length(across_align_tr)));
for iDate = 1:length(across_align_tr)
    temp = cellfun(@(x) mean(x(avg_dims)),{across_align_tr(iDate).cc});
    val_across_tr(iDate) = mean(temp);
end
% unaligned tune to tune
[val_unalign_tt] = deal(zeros(1,length(across_unalign_tt)));
for iDate = 1:length(across_unalign_tt)
    temp = cellfun(@(x) mean(x(avg_dims)),{across_unalign_tt(iDate).r});
    val_unalign_tt(iDate) = mean(temp);
end
[val_across_rr] = deal(zeros(1,length(across_align_rr)));
for iDate = 1:length(across_align_rr)
    temp = cellfun(@(x) mean(x(avg_dims)),{across_align_rr(iDate).cc});
    val_across_rr(iDate) = mean(temp);
end


% Second figure: Compare tuned to real
bins = 0:0.01:1;
subplot(1,2,1);
hold all;
[n,x] = hist(val_within_rr,bins);
h = bar(x,n/sum(n),'hist');
set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3);


% [n,x] = hist(val_across_tr,bins);
% h = bar(x,n/sum(n),'hist');
% set(h,'FaceColor','m','EdgeColor','m','FaceAlpha',0.3);


[n,x] = hist(val_across_rr,bins);
h = bar(x,n/sum(n),'hist');
set(h,'FaceColor','b','EdgeColor','b','FaceAlpha',0.3);

[n,x] = hist(val_across_tt,bins);
h = bar(x,n/sum(n),'hist');
set(h,'FaceColor','r','EdgeColor','r','FaceAlpha',0.3);


set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 1]);
xlabel('Mean of First 4 CCs')
ylabel('Count');

h = legend({'Within-RealToReal','Across-RealToReal','Across-ModelToModel','Across-ModelToModel-CID'},'Location','NorthWest');
set(h,'Box','off','FontSize',14);

% p = ranksum(val_across_rr,val_across_tt);
% title(['p= ' num2str(p,3)]);
title('CCs')


% third figure: normalized
val_norm = mean(within_vals_rr(:,avg_dims),2)';
bins = 0:0.01:1.15;

subplot(1,2,2);
hold all

[n,x] = hist(val_across_rr./val_norm,bins);
h = bar(x,n/sum(n),'hist');
set(h,'FaceColor','b','EdgeColor','b','FaceAlpha',0.3);


% [n,x] = hist(val_across_tr./val_norm,bins);
% h = bar(x,n/sum(n),'hist');
% set(h,'FaceColor','m','EdgeColor','m','FaceAlpha',0.3);

[n,x] = hist(val_across_tt./val_norm,bins);
h = bar(x,n/sum(n),'hist');
set(h,'FaceColor','r','EdgeColor','r','FaceAlpha',0.3);

set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 bins(end)]);
xlabel('Normalized CC')
ylabel('Count');

h = legend({'RealToReal','ModelToModel','ModelToModel-CID'},'Location','NorthWest');
set(h,'Box','off','FontSize',14);

% p = ranksum(val_across_rr./val_norm,val_across_tt./val_norm);
% title(['p= ' num2str(p,3)]);
title('Normalized by within-day');


if save_figs
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.fig']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.pdf']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.png']));
end


%%




