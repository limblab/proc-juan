%
%
% Compute VAF after PAs
%
%

mani_dims = 12;


all_var_ratio = [];
all_var_ratio_no_align = [];


for d = 1:length(datasets)
    
    comb_tasks = nchoosek(1:length(datasets{d}.labels),2);
    
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
        
        all_var_ratio = [all_var_ratio, var_ratio1, var_ratio2];
        
        
        % -----------------------------------------------------------------
        % Control: variance without aligning
        
        var_no_align1_mani2 = var(fr1 * w2);
        var_no_align2_mani1 = var(fr2 * w1);
        
        var_ratio1_no_align = sum(var_no_align1_mani2)/sum(vaf_orig1);
        var_ratio2_no_align = sum(var_no_align2_mani1)/sum(vaf_orig2);
        
        
        all_var_ratio_no_align = [all_var_ratio_no_align, var_ratio1_no_align, var_ratio2_no_align];
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

y_plot = ceil(max(hist_vaf_ratio))+5;
y_stats = (y_plot - max(hist_vaf_ratio))/2+max(hist_vaf_ratio);


figure,hold on
bar(x_hist(1:end-1), hist_vaf_ratio, 'histc','facecolor','k')
plot([m_vaf_ratio-sd_vaf_ratio, m_vaf_ratio+sd_vaf_ratio],[y_stats, y_stats],'k','linewidth',3)
plot(m_vaf_ratio,y_stats,'.k','markersize',32)
set(gca,'TickDir','out','FontSize',14), box off, ylim([0 y_plot]),xlim([0 100])
xlabel('Variance explained projection onto cross-task manifold')
ylabel('Percentage (%)')