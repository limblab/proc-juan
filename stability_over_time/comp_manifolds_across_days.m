%
% Control that compares the manifold orientation over days
%

function results = comp_manifolds_across_days( td, params ) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES

mani_dims           = 1:10;
min_n_units         = 15;

if nargin > 1, assignParams(who,params); end % overwrite 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STUFF

% get all pairs of sessions
sessions                = unique({td.date});
n_sessions              = length(sessions);
comb_sessions           = nchoosek(1:n_sessions,2);
n_comb_sessions         = size(comb_sessions,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THING ITSELF

princ_angles            = nan(n_comb_sessions,length(mani_dims));
vaf_ratio               = nan(n_comb_sessions,2);

for c = 1:n_comb_sessions
   
    s1                  = comb_sessions(c,1);
    s2                  = comb_sessions(c,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0. PREPARE THE DATA    
    
    % keep the channels that are common for these two sessions
    [~,td12]            = getTDidx( td, {'date',{sessions{s1},sessions{s2}} } );
    
    td12                = getCommonUnits(td12);
    
    n_units             = size(td12(1).(params.spiking_inputs{1}),2);
    
    % check that we have enough common units
    if n_units > min_n_units
    
        % split into sessions and do PCA separately
        [~,td1]      	= getTDidx( td12, {'date',sessions{s1}} );
        [~,td2]         = getTDidx( td12, {'date',sessions{s2}} );    

        [~,pca1]        = getPCA( td1, struct( ... 
                                'signals', params.spiking_inputs ) );
        [~,pca2]        = getPCA( td2, struct( ... 
                                'signals', params.spiking_inputs ) );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. PRINCIPAL ANGLES

        w1              = pca1.w(:,mani_dims);
        w2              = pca2.w(:,mani_dims);

        princ_angles(c,:) = principal_angles( w1, w2 );


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. CROSS-TASK TO WITHIN-TASK VARIANCE EXPLAINED RATIO
        % Or how much neural variance from day 1 lies within the manifold from
        % day 2, compared to the variance within the original manifold (day 1)

        sp1          	= cat(1,td1.M1_spikes);
        sp2             = cat(1,td2.M1_spikes);

        % Within-day projections
        sc1_within      = w1'*sp1';
        sc2_within      = w2'*sp2';

        within_1        = sum( var(sc1_within,1,2) ) / sum( var(sp1) );
        within_2        = sum( var(sc2_within,1,2) ) / sum( var(sp2) );

        % Cross-day ratio
        across_onto1    = sum( var(w1'*sp2',1,2) ) / sum( var(sp2) );
        across_onto2    = sum( var(w2'*sp1',1,2) ) / sum( var(sp1) );

        vaf_ratio(c,1)  = across_onto1/within_1;
        vaf_ratio(c,2)  = across_onto2/within_2;
    end
end

%%


results.vaf_ratio       = vaf_ratio;
results.princ_angles    = princ_angles;


% get number of days between sessions, to summarize the results
diff_days   = zeros(1,size(comb_sessions,1));
for c = 1:size(comb_sessions,1)
    diff_days(c) = datenum(sessions{comb_sessions(c,2)}) - datenum(sessions{comb_sessions(c,1)});
end
results.diff_days       = diff_days;
results.comb_sessions   = comb_sessions;





pas_plot                = 1:5;

dd                      = results.diff_days(~isnan(results.vaf_ratio(:,1)));
vafr                    = vaf_ratio(~isnan(results.vaf_ratio(:,1)),:);
pas                     = princ_angles(~isnan(results.vaf_ratio(:,1)),pas_plot);

m_vafr                  = nanmean(vafr,2);
m_pas                   = rad2deg(nanmean(pas(:,pas_plot),2));
sd_pas                  = rad2deg(nanstd(pas(:,pas_plot),0,2));


% linear fits
xlf                     = [0 max(dd)+1];
lf_vafr                 = polyfit( dd, m_vafr', 1 );
lf_pas                  = polyfit( dd, m_pas', 1 );
y_vafr                  = polyval(lf_vafr,xlf);
y_pas                   = polyval(lf_pas,xlf);

% plot
figure,
subplot(211), hold on
plot(xlf,y_vafr,'k','linewidth',1.5)
plot( dd, m_vafr, '.k', 'markersize', 32 )
plot( [0 max(dd)], [1 1], '-.', 'color', [.65 .65 .65],'linewidth',1.5)
ylabel('Across-day / within-day VAF')
set(gca,'TickDir','out','FontSize',14), box off
ylim([0 1.2]), xlim([0 max(dd)+1])
subplot(212), hold on
plot(xlf,y_pas,'b','linewidth',1.5)
errorbar( dd, m_pas, sd_pas, '.b', 'markersize', 32, 'linestyle', 'none' )
ylabel(['Mean top ' num2str(pas_plot(end)) ' PAs']),
ylim([0 90]),  xlim([0 max(dd)+1])
xlabel('Days between sessions')
set(gca,'TickDir','out','FontSize',14), box off
set(gcf, 'color', [1 1 1])