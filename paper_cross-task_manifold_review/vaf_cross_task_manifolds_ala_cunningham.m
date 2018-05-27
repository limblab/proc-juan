%
%
% Compute VAF after PAs
%
%

if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

mani_dims = 12;

n_shuffles = 10000; % for the control
P = 0.01; % signif threshold for the control



% preallocate
all_var_ratio = [];
all_var_ratio_no_align = [];
all_var_shuffled = [];


% do
for d = 1:length(datasets)
    
    comb_tasks = nchoosek(1:length(datasets{d}.labels),2);
    
    disp(d);
    
    % FOR EACH TASK COMBINATION
    for c = 1:size(comb_tasks,1)
        
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);
        
        
        % -----------------------------------------------------------------
        % compute PAs
        
        w1 = datasets{d}.dim_red_FR{t1}.w(:,1:mani_dims);
        w2 = datasets{d}.dim_red_FR{t2}.w(:,1:mani_dims);
        
        [tPA, U, ~, V] = principal_angles(w1,w2);
        
        
        % -----------------------------------------------------------------
        % compute VAF in original manifolds
        fr1 = datasets{d}.stdata{t1}.target{end}.neural_data.conc_smoothed_fr;
        fr2 = datasets{d}.stdata{t2}.target{end}.neural_data.conc_smoothed_fr;
        
        lv1 = fr1 * w1;
        lv2 = fr2 * w2;
        
        vaf_orig1 = var(lv1,1);
        vaf_orig2 = var(lv2,1);
        
        
        % -----------------------------------------------------------------
        % compute VAF in aligned manifolds
        
        % Project neural data of task i into aligned manifold from task j
        % to get lvi_manij
        lv1_mani2 = fr1 * w2 * V;
        lv2_mani1 = fr2 * w1 * U;
        
        vaf1_mani2 = var(lv1_mani2);
        vaf2_mani1 = var(lv2_mani1);
        
        
        % -----------------------------------------------------------------
        % Compute variance ratio
        
        var_ratio1 = sum(vaf1_mani2)/sum(vaf_orig1);
        var_ratio2 = sum(vaf2_mani1)/sum(vaf_orig2);
        
        all_var_ratio = [all_var_ratio, var_ratio1, var_ratio2]; %#ok<*AGROW>
        
        
        % -----------------------------------------------------------------
        % Variance without aligning
        
        var_no_align1_mani2 = var(fr1 * w2);
        var_no_align2_mani1 = var(fr2 * w1);
        
        var_ratio1_no_align = sum(var_no_align1_mani2)/sum(vaf_orig1);
        var_ratio2_no_align = sum(var_no_align2_mani1)/sum(vaf_orig2);
        
        
        all_var_ratio_no_align = [all_var_ratio_no_align, var_ratio1_no_align, var_ratio2_no_align];
        
        
        % -----------------------------------------------------------------
        % CONTROL: PROJECT ONTO RANDOM MANIFOLDS
        
        all_sh1 = zeros(1,n_shuffles);
        all_sh2 = zeros(1,n_shuffles);
        
        for s = 1:n_shuffles
            
            nel = numel(w1);
            nn = size(w1,1);
            
            % do for task 1
            w1_shuf = orth(randn(nn,mani_dims));
            
            var_sh_mani1 = var(fr2 * w1_shuf);
            
            var_ratio1_sh = sum(var_sh_mani1) / sum(vaf_orig2);
            
            % do for task 2
            w2_shuf = orth(randn(nn,mani_dims));
            
            var_sh_mani2 = var(fr1 * w2_shuf);
            
            var_ratio2_sh = sum(var_sh_mani2) / sum(vaf_orig1);
            
            all_sh1(s) = var_ratio1_sh;
            all_sh2(s) = var_ratio2_sh;
        end
        
        % Statistics
        th_all_sh1 = prctile(all_sh1,100-P*100);
        th_all_sh2 = prctile(all_sh2,100-P*100);
        
        all_var_shuffled = [all_var_shuffled, th_all_sh1, th_all_sh2];
    end
end


% -------------------------------------------------------------------------
% PLOT RESULTS

% convert to percentage
all_var_ratio = all_var_ratio * 100;

m_vaf_ratio = mean(all_var_ratio);
sd_vaf_ratio = std(all_var_ratio);

x_hist = 0:5:105;
hist_vaf_ratio = histcounts(all_var_ratio,x_hist)/length(all_var_ratio)*100;


% paired t-test
% [~, p] = ttest2(double(all_var_ratio),double(all_var_ratio_sh)); % this
% is unpaired
[~, pp] = ttest(double(all_var_ratio),double(all_var_ratio_sh));


% and for the control
all_var_ratio_sh = all_var_shuffled * 100;

m_vaf_sh = mean(all_var_ratio_sh);
sd_vaf_sh = std(all_var_ratio_sh);

hist_vaf_sh = histcounts(all_var_ratio_sh,x_hist)/length(all_var_ratio_sh)*100;


y_plot = max( ceil(max(hist_vaf_ratio)), ceil(max(hist_vaf_sh))) + 5;
y_stats = (y_plot - max(hist_vaf_ratio))/2+max(hist_vaf_ratio);



figure,hold on
b1 = bar(x_hist(1:end-1), hist_vaf_ratio, 'histc');
set(b1,'FaceColor','k')
b2 = bar(x_hist(1:end-1), hist_vaf_sh, 'histc');
set(b2,'FaceColor',[.7 .7 .7])
plot([m_vaf_ratio-sd_vaf_ratio, m_vaf_ratio+sd_vaf_ratio],[y_stats, y_stats],'k','linewidth',3)
plot(m_vaf_ratio,y_stats,'.k','markersize',32)
plot([m_vaf_sh-sd_vaf_sh, m_vaf_sh+sd_vaf_sh],[y_stats, y_stats],'color',[.7 .7 .7],'linewidth',3)
plot(m_vaf_sh,y_stats,'.','markersize',32,'color',[.7 .7 .7])
set(gca,'TickDir','out','FontSize',14), box off, ylim([0 y_plot]),xlim([0 100])
xlabel('Variance explained projection onto cross-task manifold (%)')
ylabel('Percentage (%)')
legend('data',['P=' num2str(P)],'Location','NorthWest'); legend boxoff
text(5,y_stats-10,['n=' num2str(length(all_var_ratio))],'Fontsize',14)
title(['unpaired t-test: P=' num2str(pp)])

