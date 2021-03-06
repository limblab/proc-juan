%% load data
clear all;
close all;
clc;

[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  false;



which_cells = 1;
avg_dims = 1:3;

monkeys = {'chewie','chewie2','jaco','mihili'};

figure('Position',[100 100 800 800]);
within_align_store = [];
for which_cells = 1:3
    [within_align, model_align, shuff_align, all_r] = deal([]);
    for iMonk = 1:length(monkeys)
        monkey = monkeys{iMonk};
        
        load([monkey '_data.mat']);
        
        pars.mani_dims = 1:3;
        pars.align_latent_params.mani_dims = 1:3;
        
        master_td = smoothSignals(master_td,struct('signals','M1_spikes','width',pars.kernel_SD,'calc_rate',true));
        master_td_all_trials = smoothSignals(master_td_all_trials,struct('signals','M1_spikes','width',pars.kernel_SD,'calc_rate',true));
        
        
        % fit tuning curves for all cells to velocity
        %   TCs are fit on whole reach, so they have the "movement information" but
        %   aren't tied to the dynamics
        dates = unique({master_td.date});
        
        for iDate = 1:length(dates)
            [~,td] = getTDidx(master_td,'date',dates{iDate});
            
            
            td = dupeAndShift(td,'vel',-3);
            
            td_temp = td;%binTD(td,'average');
            
            
            fr = getSig(td_temp,'M1_spikes');
            v = getSig(td_temp,'vel_shift');
            s = sqrt(v(:,1).^2 + v(:,2).^2);
            
            th = atan2(v(:,2),v(:,1));
            %X = [ones(size(v,1),1) sin(th) cos(th)];
            X = [ones(size(v,1),1) s.*sin(th) s.*cos(th) s];
            
            tc = zeros(size(fr,2),4);
            r = zeros(size(fr,2),1);
            for unit = 1:size(fr,2)
                [b,~,~,~,temp] = regress(fr(:,unit),X);
                r(unit) = temp(1);
                tc(unit,:) = b;
                % convert to model b0 + b1*cos(theta+b2)
                %tc(unit,:)  = [b(1); sqrt(b(2).^2 + b(3).^2); atan2(b(2),b(3))];
                
            end
            
            
            for trial = 1:length(td)
                idx = randperm(size(td(trial).M1_spikes,1));
                td(trial).vel_shuff = td(trial).vel_shift(idx,:);
                td(trial).M1_spikes_shuff = td(trial).M1_spikes(idx,:);
            end
            
            
            v = getSig(td,'vel_shift');
            s = sqrt(v(:,1).^2 + v(:,2).^2);
            
            th = atan2(v(:,2),v(:,1));
            % project kinematics through TC to get fake firing rates for all reaches
            fr_model = zeros(size(fr));
            for unit = 1:size(fr,2)
                fr_model(:,unit) = tc(unit,1) + tc(unit,4)*s + tc(unit,2).*s.*sin(th) + tc(unit,3).*s.*cos(th);
                %fr_model(:,unit) = tc(unit,1) + tc(unit,2)*cos(th + tc(unit,3));
            end
            
            
            
            v = getSig(td,'vel_shuff');
            s = sqrt(v(:,1).^2 + v(:,2).^2);
            
            th = atan2(v(:,2),v(:,1));
            % project kinematics through TC to get fake firing rates for all reaches
            fr_shuff = zeros(size(fr));
            for unit = 1:size(fr,2)
                fr_shuff(:,unit) = tc(unit,1) + tc(unit,4)*s + tc(unit,2).*s.*sin(th) + tc(unit,3).*s.*cos(th);
                %fr_model(:,unit) = tc(unit,1) + tc(unit,2)*cos(th + tc(unit,3));
            end
            
            % add modeled FR to TD
            td2 = td;
            count = 0;
            for trial = 1:length(td)
                td2(trial).M1_spikes = fr_model(count+1:count+size(td2(trial).M1_spikes,1),:);
                count = count + size(td2(trial).M1_spikes,1);
            end
            
            td3 = td;
            count = 0;
            for trial = 1:length(td)
                td3(trial).M1_spikes = fr_shuff(count+1:count+size(td3(trial).M1_spikes,1),:);
                count = count + size(td3(trial).M1_spikes,1);
            end
            
            
            % do CCA between modeled and real FRs
            
            % use only the channels that are cosine-tuned
            switch which_cells
                case 1
                    disp('subsampling neurons based on tuning');
                    neuron_idx  = r < 0.2;
                case 2
                    disp('subsampling neurons based on tuning');
                    neuron_idx  = r > 0.2;
                case 3
                    neuron_idx = true(1,size(td(1).M1_spikes,2));
            end
            
            if sum(neuron_idx) > 3
                for trial = 1:length(td)
                    td(trial).M1_spikes = td(trial).M1_spikes(:,neuron_idx);
                    td2(trial).M1_spikes = td2(trial).M1_spikes(:,neuron_idx);
                    td3(trial).M1_spikes = td3(trial).M1_spikes(:,neuron_idx);
                end
                
                td = dimReduce(td,'M1_spikes');
                td2 = dimReduce(td2,'M1_spikes');
                td3 = dimReduce(td3,'M1_spikes');
                
                model_align = [model_align, compDynamics( [td,td2], 'M1_pca', 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
                
                shuff_align = [shuff_align, compDynamics( [td,td3], 'M1_pca', 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims)];
                
                all_r = [all_r; r];
            end
            
            % subsample for within-day
            idx = getTDidx(master_td_all_trials,'date',dates{iDate});
            if sum(neuron_idx) > 3
                for trial = idx
                    master_td_all_trials(trial).M1_spikes = master_td_all_trials(trial).M1_spikes(:,neuron_idx);
                end
            else
                master_td_all_trials(idx) = [];
            end
        end
        
        
        % do within day
        within_align = [within_align, align_latent_activity_within_day(master_td_all_trials, pars.align_latent_params )];
        
    end
    
    within_align_store = [within_align_store, {within_align}];
    
    % compare results
    [real,fake, shuff] = deal(zeros(1,length(model_align)));
    for iDate = 1:length(model_align)
        temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
        real(iDate) = mean(temp);
        fake(iDate) = mean(model_align(iDate).cc(avg_dims));
        shuff(iDate) = mean(shuff_align(iDate).cc(avg_dims));
    end
    
    subplot(2,2,which_cells+1);
    hold all;
    
    bins = 0:0.02:1;
    
    [n,x] = hist(real,bins);
    h = bar(x,n,'hist');
    set(h,'FaceColor','b','EdgeColor','b','FaceAlpha',0.3);
    
    [n,x] = hist(fake,bins);
    h = bar(x,n,'hist');
    set(h,'FaceColor','r','EdgeColor','r','FaceAlpha',0.3);
    
    % [n,x] = hist(shuff,bins);
    % h = bar(x,n,'hist');
    % set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3);
    
    
    set(gca,'Box','off','TickDir','out','FontSize',14)
    set(gca,'XLim',[0 1]);
    xlabel('Mean of First 3 CCs')
    ylabel('Count');
    
    h = legend({'Within-day real','Modeled to real','Shuffled model to real'},'Location','NorthWest');
    set(h,'Box','off','FontSize',14);
    
    switch which_cells
        case 1
            title('Using only non-cosine neurons')
        case 2
            title('Using only cosine neurons')
        case 3
            title('Using all neurons')
    end
    
end

all_r(isnan(all_r)) = [];

subplot(2,2,1);
hold all;
hist(all_r,0:0.02:1);
V = axis;
plot(median(all_r)*[1 1],V(3:4),'r','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 1]);
xlabel('R^2 of cosine fit')
ylabel('Count');
title(['Median: ' num2str(median(all_r))]);



if save_figs
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_CosineTuningModel.fig']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_CosineTuningModel.pdf']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_CosineTuningModel.png']));
end




%%

% plot distribution of all neurons and only non-cosine neurons

% compare results
within_align = within_align_store{1};
[vals1] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    vals1(iDate) = mean(temp);
end
within_align = within_align_store{2};
[vals2] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    vals2(iDate) = mean(temp);
end

within_align = within_align_store{3};
[vals3] = deal(zeros(1,length(within_align)));
for iDate = 1:length(within_align)
    temp = cellfun(@(x) mean(x(avg_dims)),{within_align(iDate).aligned_info.cc});
    vals3(iDate) = mean(temp);
end


bins = 0:0.02:1;

figure;

subplot(1,2,1);
hold all;
hist(all_r,0:0.02:1);
V = axis;
plot(median(all_r)*[1 1],V(3:4),'r','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 1]);
xlabel('R^2 of cosine fit')
ylabel('Count');
title(['Median: ' num2str(median(all_r))]);


subplot(1,2,2);
hold all;

[n,x] = hist(vals3,bins);
h = bar(x,n,'hist');
set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3);

[n,x] = hist(vals1,bins);
h = bar(x,n,'hist');
set(h,'FaceColor','b','EdgeColor','b','FaceAlpha',0.3);

[n,x] = hist(vals2,bins);
h = bar(x,n,'hist');
set(h,'FaceColor','g','EdgeColor','g','FaceAlpha',0.3);

set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 1]);
xlabel('Mean of First 3 CCs')
ylabel('Count');

h = legend({'All Neurons','Untuned Neurons','Tuned Neurons'});
set(h,'Box','off','FontSize',14);



if save_figs
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.fig']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.pdf']));
    saveas(gcf,fullfile(save_dir,'Cosine Tuning Model',[ monkey '_' array '_AllVsNoncosine_Hist.png']));
end



