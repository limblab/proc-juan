%
% PA dropping channels, done properly
%

mani_dims = 12;

ds = [1 8 11 4];
tasks = [3 4 2 2];

perc_drop = [.5 .4 .3 .2 .1];
reps = 100;


% load all the data
if ~exist('datasets','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/all_manifold_datasets.mat');
end

% load significance threshold
if ~exist('angle_non_orth','var')
    load('/Users/juangallego/Documents/NeuroPlast/Data/_Dimensionality reduction/_control analyses/empirical principal angle distributions all datasets_original_submission.mat');
end


PAs_shuf = zeros(reps,length(perc_drop),mani_dims,length(ds));


for d = 1:length(ds)
   
    fr = datasets{ds(d)}.stdata{tasks(d)}.target{end}.neural_data.conc_smoothed_fr;
    
    n_units = size(fr,2);
    
    for p = 1:length(perc_drop)
        
        for r = 1:reps

            % drop a certain percentage
            rand1 = randperm(n_units);
            drop1 = rand1(1:floor(length(rand1)*perc_drop(p)));

            rand2 = randperm(n_units);
            drop2 = rand2(1:floor(length(rand2)*perc_drop(p)));

            fr1 = fr;
            fr2 = fr;

            fr1(:,drop1) = zeros(size(fr1,1),length(drop1));
            fr2(:,drop2) = zeros(size(fr2,1),length(drop2));

            w1 = pca(fr1);
            w2 = pca(fr2);

            w1 = w1(:,1:mani_dims);
            w2 = w2(:,1:mani_dims);

            PAs_shuf(r,p,:,d) = principal_angles(w1,w2);
        end
    end
end

% stats
mn_PA = mean(PAs_shuf,1);

cols = parula(length(perc_drop));


figure,
for d = 1:length(ds)
    
    ano = angle_non_orth(:,1,ds(d));
    
    subplot(1,length(ds),d), hold on
    for p = 1:length(perc_drop) 
        plot(rad2deg(squeeze(mn_PA(:,p,:,d))),'color',cols(p,:),'linewidth',1.5)
    end
    plot(ano,'color',[.65 .65 .65],'linewidth',1.5,'linestyle','-.')
    set(gcf,'color','w')
    set(gca,'FontSize',14,'TickDir','out')
    title([datasets{ds(d)}.monkey ' - ' datasets{ds(d)}.labels{tasks(d)} ' - N = ' num2str(length(datasets{ds(d)}.neural_chs))]);
    
    if d == 1
        for p = 1:length(perc_drop)
            lgnd{p} = [num2str(perc_drop(p)*100) ' % drop'];
        end
        lgnd{length(lgnd)+1} = 'P<0.001';
    end
    legend(lgnd,'Location','SouthEast'),legend boxoff
end