%
% Compute variance accounted for by neural modes projected onto the
% dimensions found with principal angles analysis
%



% what datasets to use
ds_to_use = 1:11;

% dataset type
ds_wrist = [1:3 7:9];
ds_reach = [4:6 10:11];

% manifold dimensionality
mani_dims = 12;


% define matrices to store variances: separately for each task in the pair
% --don't really know why, to compare them?
all_vars_A = [];
all_vars_B = [];
wrist_var = [];


% -------------------------------------------------------------------------
% Do!

for d = 1:length(ds_to_use)
    
    % dataset idx
    t_ds = ds_to_use(d);

    % get all tasks comparisons within this session
    comb_tasks = nchoosek(1:length(datasets{t_ds}.labels),2);

    
    % do for all tasks comparisons
    for c = 1:size(comb_tasks,1)

        % task idx within dataset
        t1 = comb_tasks(c,1);
        t2 = comb_tasks(c,2);

        % -----------------------------------------------------------------
        % 0. Get data that will be used for the analyses

        % A, B: Matrices with eigenvectors defining the manifolds
        A = datasets{t_ds}.dim_red_FR{t1}.w(:,1:mani_dims);
        B = datasets{t_ds}.dim_red_FR{t2}.w(:,1:mani_dims);

        % fra, frb: Neural activity during each task
        fra = datasets{t_ds}.stdata{t1}.target{end}.neural_data.conc_smoothed_fr;
        frb = datasets{t_ds}.stdata{t2}.target{end}.neural_data.conc_smoothed_fr;

        % -----------------------------------------------------------------
        % 1. Compute principal angles, and get SVD matrices (SVD -> U·S·V')

        % Do SVD of the inner product matrix; the diagonal of S are the cosines
        % of the principal angles
        [U, S, V] = svd(A'*B);

        % The new bases for the neural modes from task A and B are:
        Apa = A*U;
        Bpa = B*V;

        % -----------------------------------------------------------------
        % Compute the variance by projecting the neural data

        % get neural modes
        modesa = fra*Apa;
        modesb = fra*Bpa;

        % compute their normalized variance
        vara = var(modesa,1)/sum(var(modesa,1));
        varb = var(modesb,1)/sum(var(modesb,1));

        % -----------------------------------------------------------------
        % Store all variances
        
        all_vars_A = [all_vars_A; vara];
        all_vars_B = [all_vars_B; varb]; %#ok<*AGROW>
        
        if ismember(t_ds,ds_wrist)
            wrist_var = [wrist_var; 1];
        else
            wrist_var = [wrist_var; 0];
        end
    end
end


% -------------------------------------------------------------------------
% Summary plots


% combined both VAFs (all_vars_A and all_vars_B) since it doesn't really
% change

all_vars = [all_vars_A; all_vars_B];

mn_all_vars = mean(all_vars,1);
sd_all_vars = std(all_vars,0,1);

figure, hold on
plot(100*all_vars','color',[.6 .6 .6])
plot(100*mn_all_vars,'k','linewidth',2)
plot(100*(mn_all_vars+sd_all_vars),'-.k','linewidth',2)
plot(100*(mn_all_vars-sd_all_vars),'-.k','linewidth',2)
% errorbar(100*mn_all_vars,100*sd_all_vars,'.k','linestyle','none','linewidth',1.5,'markersize',30)
set(gca,'TickDir','out','FontSize',14), box off
xlabel('Principal-angle projected neural mode')
ylabel('Neural variance explained (%)')


% % separating by tasks
% figure, hold on
% plot(100*all_vars(wrist_var==1,:)','color',[.6 .6 .6])
% plot(100*all_vars(wrist_var==0,:)','color',[1 .6 0])
% % errorbar(100*mn_all_vars,100*sd_all_vars,'.k','linestyle','none','linewidth',1.5,'markersize',30)
% set(gca,'TickDir','out','FontSize',14), box off
% xlabel('Neural mode (after computing principal angles)')
% ylabel('Neural variance explained (%)')
