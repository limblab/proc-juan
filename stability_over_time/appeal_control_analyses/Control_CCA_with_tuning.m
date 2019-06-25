%% load data
clear all;
close all;
clc;


[save_dir, ~] = get_computer_paths();
save_dir = fullfile(save_dir,'_Appeal_Results');

save_figs =  false;

min_fr = 0;
r2_cutoff = 0.6;
do_speed = false;

mani_dims = 1:8;
avg_dims = 1:4;

monkeys = {'chewie','chewie2','jaco','mihili'};
array =  'M1';

[monk_labels, within_align, all_to_tune, all_to_untune, tune_to_untune, all_r,decode_results] = deal([]);
for iMonk = 1:length(monkeys)
    monkey = monkeys{iMonk};
    
    load([monkey '_data.mat']);
    
    pars.mani_dims = mani_dims;
    pars.align_latent_params.mani_dims = mani_dims;
    
    master_td_all_trials = smoothSignals(master_td_all_trials,struct('signals','M1_spikes','width',pars.kernel_SD,'calc_rate',true));
    
    % fit tuning curves for all cells to velocity
    %   TCs are fit on whole reach, so they have the "movement information" but
    %   aren't tied to the dynamics
    dates = unique({master_td_all_trials.date});
    
    within_td = [];
    for iDate = 1:length(dates)
        [~,td] = getTDidx(master_td_all_trials,'date',dates{iDate});
        td = removeBadNeurons(td,struct('min_fr',min_fr));
        
        if do_speed
            td_temp = td;
        else
            td_temp = binTD(td,'average');
        end
        
        fr = getSig(td_temp,'M1_spikes');
        v = getSig(td_temp,'vel');
        s = sqrt(v(:,1).^2 + v(:,2).^2);
        
        th = atan2(v(:,2),v(:,1));
        
        if ~do_speed % no speed
            disp('Tuning with no speed');
            X = [ones(size(v,1),1) sin(th) cos(th)];
            
            tc = zeros(size(fr,2),3);
            r = zeros(size(fr,2),1);
            for unit = 1:size(fr,2)
                [b,~,~,~,temp] = regress(fr(:,unit),X);
                r(unit) = temp(1);
                tc(unit,:) = b;
                % convert to model b0 + b1*cos(theta+b2)
                %tc(unit,:)  = [b(1); sqrt(b(2).^2 + b(3).^2); atan2(b(2),b(3))];
                
            end
        else
            disp('Tuning with speed');
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
        end
        all_r = [all_r; r];
        
        % do CCA between tuning conditions
        % tuned
        neuron_idx  = r > r2_cutoff;
        if sum(neuron_idx) > mani_dims(end) && sum(~neuron_idx) > mani_dims(end)
            within_td = [within_td, td];
            
            % set up our extra tds
            td_temp=[];
            count  = 0;
            while length(td_temp) < 120
                count = count + 1;
                trial_idx = randperm(length(td));
                td_all = td(trial_idx(floor(length(trial_idx)/2)+1:end));
                td_tune = td(trial_idx(1:floor(length(trial_idx)/2)));
                td_untune = td(trial_idx(1:floor(length(trial_idx)/2)));
                
                % give fake dates to the extras
                for trial = 1:length(td_tune)
                    td_tune(trial).date = '01-01-2017';
                    td_untune(trial).date = '02-01-2017';
                end
                
                td_temp = equalNbrTrialsSessions([td_all, td_tune, td_untune]);
            end
            
            fuck = unique({td_temp.date});
            [~,td_tune]  = getTDidx(td_temp,'date',fuck{1});
            [~,td_untune]  = getTDidx(td_temp,'date',fuck{2});
            [~,td_all]  = getTDidx(td_temp,'date',fuck{3});
            
            
            
            for trial = 1:length(td_tune)
                td_tune(trial).M1_spikes = td_tune(trial).M1_spikes(:,neuron_idx);
            end
            %  untuned
            neuron_idx  = r <= r2_cutoff;
            for trial = 1:length(td_untune)
                td_untune(trial).M1_spikes = td_untune(trial).M1_spikes(:,neuron_idx);
            end
            
            td_all = dimReduce(td_all,'M1_spikes');
            td_tune = dimReduce(td_tune,'M1_spikes');
            td_untune = dimReduce(td_untune,'M1_spikes');
            
            all_to_tune = [all_to_tune, compDynamics( [td_all,td_tune], 'M1_pca', 1:length(td_all), length(td_all)+1:length(td_all)+length(td_tune),pars.mani_dims )];
            all_to_untune = [all_to_untune, compDynamics( [td_all,td_untune], 'M1_pca', 1:length(td_all), length(td_all)+1:length(td_all)+length(td_untune),pars.mani_dims)];
            tune_to_untune = [tune_to_untune, compDynamics( [td_tune,td_untune], 'M1_pca', 1:length(td_tune), length(td_tune)+1:length(td_tune)+length(td_untune),pars.mani_dims)];
            monk_labels = [monk_labels, {monkey}];
        end
        
        
    end
    
    within_align = [within_align, align_latent_activity_within_day(within_td,pars.align_latent_params)];
end


%% compare results

[vals0,vals1,vals2, vals3] = deal(zeros(1,length(all_to_tune)));
for iDate = 1:length(all_to_tune)
    vals0(iDate) = mean(within_align(iDate).cc_m(avg_dims));
    vals1(iDate) = mean(all_to_tune(iDate).cc(avg_dims));
    vals2(iDate) = mean(all_to_untune(iDate).cc(avg_dims));
        vals3(iDate) = mean(tune_to_untune(iDate).cc(avg_dims));
end

% vals1 = vals1./vals0;
% vals2 = vals2./vals0;



bins = 0:0.05:1;

figure;
all_r(isnan(all_r)) = [];
subplot(1,2,1);
hold all;
hist(all_r,0:0.025:1);
V = axis;
plot(r2_cutoff*[1 1],V(3:4),'r','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14)
set(gca,'XLim',[0 bins(end)]);
xlabel('R^2 of cosine fit')
ylabel('Count');
title(['Cutoff: ' num2str(r2_cutoff)]);


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
set(gca,'XLim',[0 bins(end)]);
xlabel('Mean of First 4 CCs')
ylabel('Count');

h = legend({'All to Tuned Neurons','All to Untuned Neurons','Tuned to Untuned Neurons'});
set(h,'Box','off','FontSize',14);

if save_figs
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned.fig']));
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned.pdf']));
    saveas(gcf,fullfile(save_dir,'CCA for Untuned',[ 'AllMonkeys_' array '_CCAforUntuned.png']));
end




% %
% % % %%
% % % c = cov(getSig(td,'M1_spikes'));
% % % % c = log10(abs(c));
% % % c = zscore(c);
% % % % c(eye(size(c)) == 1) = 0;
% % %
% % % figure;
% % % hold all;
% % % % now sort by  tuning
% % % [~,I]  = sort(r);
% % % imagesc(c(I,I));
% % % axis square;
% % % axis tight;
% % % set(gca,'Box','off','TickDir','out','FontSize',14);
% % %
% % % % colorbar;
% %
% % %%
% % % compute decoding results
% % dec_within = cellfun(@mean,{decode_results.withinR2_m});
% % dec_across = cellfun(@mean,{decode_results.acrossR2});
% % dec_ctrl = cellfun(@mean,{decode_results.ctrlR2});
% %
% % % plot decoding
% % figure
% % hold all;
% % boxplot([dec_ctrl; dec_within; dec_across]');
% % set(gca,'Box','off','TickDir','out','FontSize',14);
% % set(gca,'XTick',1:3,'XTickLabel',{'Untuned','Tuned','Across'});
% % set(gca,'YLim',[0 1]);
% % ylabel('R^2');
% % title('Decoder Perf.');




