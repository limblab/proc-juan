%
% PLOTS PREDICTIONS!!!
%
% Needs to be run pausing in decode_across_Days. Fixing this is easy
% enough, but we don't have time!
%


% COLORS
col_xpred = [139 0 139]/255;
col_ypred = [255,20,147]/255;
col_xpredwithin = [255 0 0]/255;
col_ypredwithin = [255 165 0]/255;
col_xspikes = [0 100 0]/255;
col_yspikes = [173 255 47]/255;
col_x = [0 0 0];
col_y = [.6 .6 .6];


if ~exist('Yhat_across','var')
    error('You need to run this within decode_across days... or save the predictions');
end



% SINGLE TRIAL PREDICTIONS


% WHAT WE HAVE:
% td2 has the predictions with the model trained on day 1 (in
% linmodel_default)
% mod_info has the model

% some definitions
binsize = td2(1).bin_size;
trialdur = size(td2(1).pos,1);
trialsep = binsize*2;


% Since we also want to compute the within-day predictions and the
% cross-days predictions based on spikes, do that!
[~,this_td1] = getTDidx(td, 'date', sessions{comb_sessions(c,1)} );
[~,this_td2] = getTDidx(td, 'date', sessions{comb_sessions(c,2)} );

% In the main function we shuffle the trials, apply it here
this_td1 = this_td1(idx_rnd);
this_td2 = this_td2(idx_rnd);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Compute the within day model based on PCs, as a ceiling --decode_across_day
% only does it crossvalidated
within_day_mod_params = mod_params;
within_day_mod_params.in_signals{1} = [mod_params.in_signals{1}(1:end-12) '_shift'];

n_pca_dims = size(this_td2(1).(within_day_mod_params.in_signals{1}(1:end-6)),2);

% dupe and shift
this_td2_within = dupeAndShift(this_td2,{within_day_mod_params.in_signals{1}(1:end-6),hist_bins});

shifted_dims_use = [];
for hb = 1:hist_bins
    shifted_dims_use = [shifted_dims_use, mani_dims+(hb-1)*n_pca_dims];
end
within_day_mod_params.in_signals{2} = shifted_dims_use;

% train model and compute R2 (this is not cross-validated!)
[this_td2_within, within_day_mod_info] = getModel(this_td2_within,within_day_mod_params);
[R2_within_notxval, ~] = testModel(this_td2_within,within_day_mod_info);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Test a decoder based on spikes

% optional: "unsort" channels

if unsort_chs_to_pred
    
    this_td1_unsorted = stripSpikeSorting( this_td1 );
    this_td2_unsorted = stripSpikeSorting( this_td2 );
end

% a) Get the channels that have threshold crossings in both datasets
both_tds = [this_td1_unsorted, this_td2_unsorted];
both_tds = getCommonUnits(both_tds);

% b) Create TD structs that only have common units
[~,this_td1_unsorted] = getTDidx(both_tds, 'date', sessions{comb_sessions(c,1)} );
[~,this_td2_unsorted] = getTDidx(both_tds, 'date', sessions{comb_sessions(c,2)} );

% c) Set model params and inputs
spike_in    = [manifold(1:end-4) '_spikes'];
n_units     = size(td1(1).(spike_in),2);
spike_mod_params = mod_params;
spike_mod_params.in_signals = {[spike_in '_shift'],1:1:n_units*hist_bins};

% d) Duplicate and shift to have history
this_td1_unsorted = dupeAndShift(this_td1_unsorted,{spike_in,hist_bins});
this_td2_unsorted = dupeAndShift(this_td2_unsorted,{spike_in,hist_bins});

% e) compute model and test it on the same session
[this_td1_unsorted, spike_mod_info] = getModel(this_td1_unsorted,spike_mod_params);
[R2_spikes,~] = testModel(this_td1_unsorted,spike_mod_info);

% TEST THE SPIKE MODEL on the second session
[R2_spikes_across, this_td2_unsorted] = testModel(this_td2_unsorted,spike_mod_info);


% -------------------------------------------------------------------------
figure,
subplot(211),hold on; subplot(212), hold on;
for tri = 1:length(td2)
    
    tm = [0:binsize:(trialdur-1)*binsize] + (tri-1)*(trialdur)*binsize;
    if tri > 1, tm = tm + trialsep; end
    
    subplot(2,1,1),hold on
    plot(tm,td2(tri).(mod_info.out_signals{1})(:,1),'color',col_x,'linewidth',1)
    plot(tm,td2(tri).linmodel_default(:,1),'color',col_xpred,'LineWidth',1)
    plot(tm,this_td2_within(tri).linmodel_default(:,1),'color',col_xpredwithin,'LineWidth',1)
    plot(tm,this_td2_unsorted(tri).linmodel_default(:,1),'color',col_xspikes,'LineWidth',1)
    
    subplot(2,1,2),hold on
    plot(tm,td2(tri).(mod_info.out_signals{1})(:,2),'color',col_y,'linewidth',1)
    plot(tm,td2(tri).linmodel_default(:,2),'color',col_ypred,'LineWidth',1)
    plot(tm,this_td2_within(tri).linmodel_default(:,2),'color',col_ypredwithin,'LineWidth',1)
    plot(tm,this_td2_unsorted(tri).linmodel_default(:,2),'color',col_yspikes,'LineWidth',1)
    
    if tri == length(td2)
        subplot(211),set(gca,'FontSize',14,'TickDir','out'); box off;
        subplot(212),set(gca,'FontSize',14,'TickDir','out'); box off;
        set(gcf, 'color', [1 1 1])
        subplot(211),ylabel(['X ' mod_params.out_signals])
        subplot(212),ylabel(['Y ' mod_params.out_signals])
        yl = ylim;
        % subplot(211),text(tm(end-length(tm)+3),yl(2)*3/4,['R^2_{across} = ' num2str(R2diff(1),2)])
        xlabel('Time (s)')
        % subplot(212),text(tm(end-length(tm)+3),yl(2)*3/4,['R^2_{across} = ' num2str(R2diff(2),2)])
        for sp = 1:2
            subplot(2,1,sp)
            legend('Actual',['Across-day aligned R^2=' num2str(R2diff(sp),sp)],...
                ['Within-day R^2=' num2str(R2_within_notxval(sp),sp)],['Spikes R^2=' num2str(R2_spikes_across(sp),sp)],...
                'Location','NorthWest'); legend boxoff
        end
        subplot(211)
        title([td(1).monkey ' ' sessions{comb_sessions(c,1)} ' vs ' sessions{comb_sessions(c,2)}])
    end
end

