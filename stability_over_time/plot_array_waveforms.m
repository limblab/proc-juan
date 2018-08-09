clear;
close all;
clc;

which_date = 2;

date0 = '2016-09-09';
date1 = '2016-10-05';
date2 = '2016-10-21';

switch which_date
    case 1
        load(['/Users/mattperich/Desktop/Chewie/CDS/' date1 '/Chewie_CO_FF_BL_10052016_001.mat']);
    case 2
        load(['/Users/mattperich/Desktop/Chewie/CDS/' date2 '/Chewie_CO_CS_BL_10212016_001.mat'])
end

d_day(1) = datenum(date1) - datenum(date0);
d_day(2) = datenum(date2) - datenum(date0);


%%
close all;

% Loop along units and plot mean of all waveforms in the appropriate location
unit_colors = {'k','b','r','g','c','m','y'};
figure; subplot1(10,10);

array_idx = find(strcmpi({cds.units.array},'M1'));
chans = unique([cds.units(array_idx).chan]);

for i = 1:length(chans)
    unit_idx = find([cds.units.chan] == chans(i) & strcmpi({cds.units.array},'M1'));
    
    rows = [cds.units(unit_idx).rowNum];
    cols = [cds.units(unit_idx).colNum];
    
    row = unique(rows);
    col = unique(cols);
    
    if length(row) > 1 || length(col) > 1
        error('not same');
    end
    
    did_anything = false;
    for j = unit_idx
        if cds.units(j).ID > 0 && cds.units(j).ID ~= 255            
            subplot1(10*(row-1)+col); hold all;
            plot(mean(cds.units(j).spikes.wave,1),'-','LineWidth',1,'Color',unit_colors{cds.units(j).ID});
            did_anything = true;
        end
    end
    if did_anything
        V = axis;
        text(0.1*V(2),0.9*V(4),['elec ' num2str(chans(i))]);
    end
end

for i = 1:10
    for j = 1:10
        subplot1(10*(i-1)+j);
        set(gca,'Box','off','TickDir','out','XTick',[],'YTick',[]);
        ax = gca; 
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
        axis('square');
        
        if 10*(i-1)+j == 1
            title(['Day ' num2str(d_day(which_date))]);
        end
    end
end


% % clear;
% % close all;
% % clc;
% % 
% % [save_dir, data_dir] = get_computer_paths();
% % save_them = false;
% % 
% % pars.monkey         = 'Chewie'; % 'chewie2'; 'chewie'; 'mihili'; 'han'; 'chips'; 'jaco'
% % pars.array          = 'M1';
% % pars.spiking_inputs = {'M1_spikes'}; % {'PMd_spikes'}; {'M1_spikes'}; {'S1_spikes'}
% % 
% % date0 = '2016-09-09';
% % date1 = '2016-10-05';
% % date2 = '2016-10-21';
% % 
% % %% --------------------------------------------------------------------------
% % % load CDS for first day
% % load(['/Users/mattperich/Desktop/Chewie/CDS/' date1 '/Chewie_CO_FF_BL_10052016_001.mat']);
% % cds1 = cds;
% % % load CDS for second day
% % load(['/Users/mattperich/Desktop/Chewie/CDS/' date2 '/Chewie_CO_CS_BL_10212016_001.mat']);
% % cds2 = cds;
% % clear cds;
% % 
% % d_day1 = datenum(date1) - datenum(date0);
% % d_day2 = datenum(date2) - datenum(date0);
% % 
% % 
% % %% --------------------------------------------------------------------------
% % % get the waveforms for the electrode
% % idx1 = find(td1(1).M1_unit_guide(:,1) == elec);
% % idx2 = find(td2(1).M1_unit_guide(:,1) == elec);
% % 
% % wf1 = cell(1,length(idx1));
% % for i = 1:length(idx1)
% %     idx = [cds1.units.chan] == td1(1).M1_unit_guide(idx1(i),1) & ...
% %         [cds1.units.ID] == td1(1).M1_unit_guide(idx1(i),2) & ...
% %         strcmpi({cds1.units.array},'M1');
% %     wf1{i} = cat(1,cds1.units(idx).spikes.wave);
% % end
% % 
% % wf2 = cell(1,length(idx2));
% % for i = 1:length(idx2)
% %     idx = [cds2.units.chan] == td2(1).M1_unit_guide(idx2(i),1) & ...
% %         [cds2.units.ID] == td2(1).M1_unit_guide(idx2(i),2) & ...
% %         strcmpi({cds2.units.array},'M1');
% %     wf2{i} = cat(1,cds2.units(idx).spikes.wave);
% % end
% % 
