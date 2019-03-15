%% load data
clear all;
close all;
clc;

mani_dims = 1:10;
avg_dims = 1:4;

monkeys = {'chewie','chewie2','jaco','mihili'};

r_vals = 0.05:0.05:0.3;

res_out = cell(length(r_vals),3);
for val = 1:length(r_vals)
    
    [within_align, all_to_tune, all_to_untune, tune_to_untune, all_r,decode_results] = deal([]);
    for iMonk = 1:length(monkeys)
        monkey = monkeys{iMonk};
        
        load([monkey '_data.mat']);
        
        pars.mani_dims = mani_dims;
        pars.align_latent_params.mani_dims = mani_dims;
        
        master_td = smoothSignals(master_td,struct('signals','M1_spikes','width',pars.kernel_SD,'calc_rate',true));
        master_td_all_trials = smoothSignals(master_td_all_trials,struct('signals','M1_spikes','width',pars.kernel_SD,'calc_rate',true));
        
        
        % fit tuning curves for all cells to velocity
        %   TCs are fit on whole reach, so they have the "movement information" but
        %   aren't tied to the dynamics
        dates = unique({master_td.date});
        
        for iDate = 1:length(dates)
            [~,td] = getTDidx(master_td,'date',dates{iDate});
            
            td_temp = td;%binTD(td,'average');
            
            fr = getSig(td_temp,'M1_spikes');
            v = getSig(td_temp,'vel');
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
            all_r = [all_r; r];
            
            % do CCA between tuning conditions
            td2 = td;
            td3 = td;
            neuron_idx = true(1,size(td(1).M1_spikes,2));
            for trial = 1:length(td)
                td(trial).M1_spikes = td(trial).M1_spikes(:,neuron_idx);
            end
            neuron_idx  = r >= r_vals(val);
            if sum(neuron_idx) > 10 && sum(~neuron_idx) > 10
                for trial = 1:length(td2)
                    td2(trial).M1_spikes = td2(trial).M1_spikes(:,neuron_idx);
                end
                neuron_idx  = r < r_vals(val);
                for trial = 1:length(td3)
                    td3(trial).M1_spikes = td3(trial).M1_spikes(:,neuron_idx);
                end
                
                [td,pca_info] = dimReduce(td,'M1_spikes');
                td2 = dimReduce(td2,'M1_spikes');
                td3 = dimReduce(td3,'M1_spikes');
                
                all_to_tune = [all_to_tune, compDynamics( [td,td2], 'M1_pca', 1:length(td), length(td)+1:length(td)+length(td2),pars.mani_dims )];
                all_to_untune = [all_to_untune, compDynamics( [td,td3], 'M1_pca', 1:length(td), length(td)+1:length(td)+length(td3),pars.mani_dims)];
                tune_to_untune = [tune_to_untune, compDynamics( [td2,td3], 'M1_pca', 1:length(td2), length(td2)+1:length(td2)+length(td3),pars.mani_dims)];
                
            end
            
        end
        
    end
    
    [vals1,vals2, vals3] = deal(zeros(1,length(all_to_tune)));
    for iDate = 1:length(all_to_tune)
        vals1(iDate) = mean(all_to_tune(iDate).cc(avg_dims));
        vals2(iDate) = mean(all_to_untune(iDate).cc(avg_dims));
        vals3(iDate) = mean(tune_to_untune(iDate).cc(avg_dims));
    end
    
    
    res_out{val,1} = vals1;
    res_out{val,2} = vals2;
    res_out{val,3} = vals3;
    
end


%%

figure;
hold all;
m = cellfun(@mean,res_out);
s = cellfun(@std,res_out);%./sqrt(cellfun(@length,res_out));
plot(r_vals,m,'LineWidth',3);
plot_colors = parula;
for i = 1:3
    h = patch([r_vals, fliplr(r_vals)], [m(:,i)-s(:,i); flipud(m(:,i)+s(:,i))]',plot_colors(i,:));
    set(h,'FaceColor',plot_colors(i,:),'EdgeColor',plot_colors(i,:),'FaceAlpha',0.5);
end

set(gca,'Box','off','TickDir','out','FontSize',14);


%% compare results

[vals1,vals2, vals3] = deal(zeros(1,length(all_to_tune)));
for iDate = 1:length(all_to_tune)
    vals1(iDate) = mean(all_to_tune(iDate).cc(avg_dims));
    vals2(iDate) = mean(all_to_untune(iDate).cc(avg_dims));
    vals3(iDate) = mean(tune_to_untune(iDate).cc(avg_dims));
end



bins = 0:0.02:1;

figure;
all_r(isnan(all_r)) = [];
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

[n,x] = hist(vals1,bins);
h = bar(x,n,'hist');
set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3);

[n,x] = hist(vals2,bins);
h = bar(x,n,'hist');
set(h,'FaceColor','b','EdgeColor','b','FaceAlpha',0.3);

% [n,x] = hist(vals3,bins);
% h = bar(x,n,'hist');
% set(h,'FaceColor','g','EdgeColor','g','FaceAlpha',0.3);

set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 1]);
xlabel('Mean of First 4 CCs')
ylabel('Count');

h = legend({'All to Tuned Neurons','All to Untuned Neurons','Tuned to Untuned Neurons'});
set(h,'Box','off','FontSize',14);


