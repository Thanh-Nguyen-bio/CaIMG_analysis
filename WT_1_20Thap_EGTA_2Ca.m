%% Visualize time series data from table
% AI generated Load Excel data
% Visualization by Mike X Cohen, sincxpress.com (source code)
% Slightly modification in normalization step by Thanh
close all, clear, clc
%% Load Excel time series data with multiple ROIs
filename = 'WT-1.csv';

% Read the table
dataTbl = readtable(filename);

% Display column names
disp('Available column names:');
disp(dataTbl.Properties.VariableNames')

% Identify the time column
timeColIdx = find(contains(dataTbl.Properties.VariableNames, 'Time', 'IgnoreCase', true), 1);

if isempty(timeColIdx)
    error('No time column found!');
end

% Extract time vector
timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);
timevec = timevec * 60;
npnts = length(timevec);

% Estimate sampling rate

marks= find(abs(diff(round(diff(timevec))))>=1.5) +2;

timevec = dataTbl{:, timeColIdx};
timevec = rmmissing(timevec);

% Identify fluorescence columns (exclude time column)
fluorColsIdx = find(contains(dataTbl.Properties.VariableNames, 'Ch1', 'IgnoreCase', true));

% Remove time column if accidentally matched by "Avg"
fluorColsIdx(fluorColsIdx == timeColIdx) = [];

if isempty(fluorColsIdx)
    error('No fluorescence data columns found!');
end

% Initialize fluorescence data matrix: [numROIs x numTimePoints]
neuronts = dataTbl{:, fluorColsIdx}';
neuronts = rmmissing(neuronts,2);
nROIs = size(neuronts, 1);

%% define events
evnt= {'TG 20 \muM'} ;

%% (Option) Modify timevec
timevec(1:60) = timevec(1:60)*5
timevec(61:510) = (timevec(61:510)-60)*2+300
timevec(511:end) = (timevec(511:end)-510)*5+1200

%% Visualization

% Plot an example trace
figure(5), clf
subplot(511)
plot(timevec, neuronts(9,:))
ylabel('Brightness (a.u.)')
set(gca, 'xlim', timevec([50 end]), 'xticklabel', [])
box off

% Show heatmap
subplot(5,1,2:5)
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),neuronts,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca, 'xlim', timevec([50 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar

% Plot all traces
figure(6), clf
plot(timevec, neuronts)
ylabel('Brightness (a.u.)')
title('Fluorescence of all ROIs')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([50 end]))

%% Backup raw data before normalization
neuronts_raw = neuronts;

%% get reference value of Tomato
reffluorColsIdx = find(contains(dataTbl.Properties.VariableNames, 'Ch2', 'IgnoreCase', true));
% Identify Tomato intensity columns (exclude time column)
refneuronts = dataTbl {:,reffluorColsIdx}';
refneuronts = rmmissing(refneuronts,2);
%% Visualization Tomato


% Show heatmap
figure(7), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),refneuronts,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca, 'xlim', timevec([50 end]))
view(2)
xlabel('Time (min.)')
ylabel('ROI number')
colorbar

% Plot all traces
figure(8), clf
plot(timevec, refneuronts)
ylabel('Brightness (a.u.)')
title('Tomato fluorescence of all ROIs')
xlabel('Time (min.)')
set(gca, 'xlim', timevec([50 end]))



%% Normalized to tomato 
neuronts = bsxfun(@rdivide,neuronts,refneuronts);

%%
figure(9), clf
title('Normalized to Tomato')
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),neuronts,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca,'xlim', timevec([1 end]))
xlabel('Time (min.)')
ylabel('ROI number')
colorbar
view(2)
figure(10), clf
plot(timevec, neuronts)
ylabel('Normalized to tdTomato')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)')
set(gca, 'xlim',timevec([1 end]))
disp('Select 3 time point of media changing and press ENTER to continue')
[new_marks_x,~] = ginput(3);


pause



%% Convert to dF/F

neuronts = bsxfun(@rdivide, neuronts, mean(neuronts(:,marks(1)-12:marks(1)), 2))-1;



%% dF/F visualization
figure(11), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),neuronts,'LineStyle','none')
set(gca,'ylim',[-1.5 nROIs+2.5],'FontSize',16)
set(gca,'ytick',(1:5:nROIs))
set(gca,'xlim', timevec([1 end]),'FontSize',16)
colorbar
set(gca,'clim',[-1 10],'fontsize',16)
set(gca,'YDir','reverse')
%xline(timevec(marks),'LineWidth',1.5)
grid off
view(2)

%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(end)],[-0.1 -0.1],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),-0.8,evnt, ...
    'HorizontalAlignment','center','FontSize',12)

%drawing 3 box and text to describe calcium concentration
rectangle('Position',[new_marks_x(1) nROIs+0.5 new_marks_x(2)-new_marks_x(1) 2],'FaceColor','#FFFF00')
text(new_marks_x(1)+0.5*(new_marks_x(2)-new_marks_x(1)),nROIs+1.5, ...
    '$${\ 2 mM \ Ca^{2+}}$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(2) nROIs+0.5 new_marks_x(3)-new_marks_x(2) 2],'FaceColor','#E0E0E0')
text((new_marks_x(2)+0.5*(new_marks_x(3)-new_marks_x(2))),nROIs+1.5, ...
    '$${\ Ca^{2+} \ free }$$ ','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(3) nROIs+0.5 timevec(end)-new_marks_x(3) 2],...
    'FaceColor','#FFFF00 ')
text(new_marks_x(3)+0.5*(timevec(end)-new_marks_x(3)),nROIs+1.5, ...
    '$${\ 2 mM \ Ca^{2+} }$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)


figure(12), clf
plot(timevec, neuronts)
ylabel('\DeltaF/F_0')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)')
set(gca,'xlim', timevec([1 end]),'FontSize',16)
set(gca, 'ylim', [-1.5 17])
set(gca,'YTick',(0:1:17))
set(gca,'TickDir','none')
%xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold','Interpreter','latex')
box off
%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(end)],[15 15],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),15.5,evnt, ...
    'HorizontalAlignment','center','FontSize',12)

%drawing 3 box and text to describe calcium concentration
rectangle('Position',[new_marks_x(1) -1.5 new_marks_x(2)-new_marks_x(1) 1],'FaceColor','#FFFF00')
text(new_marks_x(1)+0.5*(new_marks_x(2)-new_marks_x(1)),-1, ...
    '$${\ 2 mM \ Ca^{2+}}$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(2) -1.5 new_marks_x(3)-new_marks_x(2) 1],'FaceColor','#E0E0E0')
text((new_marks_x(2)+0.5*(new_marks_x(3)-new_marks_x(2))),-1, ...
    '$${\ Ca^{2+} \ free }$$ ','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(3) -1.5 timevec(end)-new_marks_x(3) 1],...
    'FaceColor','#FFFF00 ')
text(new_marks_x(3)+0.5*(timevec(end)-new_marks_x(3)),-1, ...
    '$${\ 2 mM \ Ca^{2+} }$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)





%% (Optional) all cell time tracing
for i=1:nROIs
    figure(100+i), clf
    plot(timevec, neuronts(i,:),'-k','LineWidth',0.25)
    ylabel('\DeltaF/F_0')
    xlabel('Time (min.)')
    title(sprintf('Fluorescence of ROIs %d (normalized)', i))
    set(gca, 'xlim', timevec([marks(1) end]))
    set(gca, 'ylim', [-1.5 10])
    set(gca,'YTick',(0:1:10))
    %xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold')
   
    %drawing 1 lines to indicate drugs treatment on top of graph
    line([timevec(marks(1)) timevec(end)],[15 15],'color','k','Linewidth',2)
    text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),15.5,evnt, ...
        'HorizontalAlignment','center','FontSize',12)

      
    %drawing 3 box and text to describe calcium concentration
    rectangle('Position',[new_marks_x(1) -1.5 new_marks_x(2)-new_marks_x(1) 1],'FaceColor','#FFFF00')
    text(new_marks_x(1)+0.5*(new_marks_x(2)-new_marks_x(1)),-1, ...
        '$${\ 2 mM \ Ca^{2+}}$$','Interpreter','latex', ...
        'HorizontalAlignment','center','FontSize',12)
    
    rectangle('Position',[new_marks_x(2) -1.5 new_marks_x(3)-new_marks_x(2) 1],'FaceColor','#E0E0E0')
    text((new_marks_x(2)+0.5*(new_marks_x(3)-new_marks_x(2))),-1, ...
        '$${\ Ca^{2+} \ free }$$ ','Interpreter','latex', ...
        'HorizontalAlignment','center','FontSize',12)
    
    rectangle('Position',[new_marks_x(3) -1.5 timevec(end)-new_marks_x(3) 1],...
        'FaceColor','#FFFF00 ')
    text(new_marks_x(3)+0.5*(new_marks_x(3)-new_marks_x(2)),-1, ...
        '$${\ 2 mM \ Ca^{2+} }$$','Interpreter','latex', ...
        'HorizontalAlignment','center','FontSize',12)

end

%% Ensure timevec is a column vector
% Force all vectors to be columns
timevec = timevec(:);
meanRaw = mean(neuronts_raw, 1)';
semRaw  = std(neuronts_raw, 0, 1)' / sqrt(size(neuronts_raw, 1));
meanNorm = mean(neuronts, 1)';
semNorm  = std(neuronts, 0, 1)' / sqrt(size(neuronts, 1));

% Plot both traces with shaded error regions
figure(13), clf

% Raw
subplot(7,1,2:4)
hold on
fill([timevec; flipud(timevec)], ...
     [meanRaw+semRaw; flipud(meanRaw-semRaw)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(timevec, meanRaw, 'b', 'LineWidth', 1.5)
ylabel('Raw Brightness (a.u.)')
title('Raw Trace with SEM')
set(gca, 'xlim', timevec([1 end]),'fontsize',12)
set(gca,'xtick',[])
set(gca,'ylim',[0 max(meanRaw)+13])
%xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold')

%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(end)],[max(meanRaw)+5 max(meanRaw)+5],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),max(meanRaw)+7,evnt, ...
    'HorizontalAlignment','center','FontSize',12)
box off

% Normalized
subplot(7,1,5:7)
hold on
fill([timevec; flipud(timevec)], ...
     [meanNorm+semNorm; flipud(meanNorm-semNorm)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(timevec, meanNorm, 'r', 'LineWidth', 1.5)
ylabel('\DeltaF/F_0')
xlabel('Time (min)')
title('\DeltaF/F_0 Trace with SEM')
set(gca, 'xlim', timevec([1 end]) , ...
    'ylim', [-1.5 10], ...
    'YTick',(0:1:10), ...
    'fontsize',12)

%xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold')



%drawing 3 box and text to describe calcium concentration
rectangle('Position',[new_marks_x(1) -1.5 new_marks_x(2)-new_marks_x(1) 1],'FaceColor','#FFFF00')
text(new_marks_x(1)+0.5*(new_marks_x(2)-new_marks_x(1)),-1, ...
    '$${\ 2 mM \ Ca^{2+}}$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(2) -1.5 new_marks_x(3)-new_marks_x(2) 1],'FaceColor','#E0E0E0')
text((new_marks_x(2)+0.5*(new_marks_x(3)-new_marks_x(2))),-1, ...
    '$${\ Ca^{2+} \ free }$$ ','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)

rectangle('Position',[new_marks_x(3) -1.5 timevec(end)-new_marks_x(3) 1],...
    'FaceColor','#FFFF00 ')
text(new_marks_x(3)+0.5*(timevec(end)-new_marks_x(3)),-1, ...
    '$${\ 2 mM \ Ca^{2+} }$$','Interpreter','latex', ...
    'HorizontalAlignment','center','FontSize',12)


box off


%% calculate thapsigargin's AUC 
%define thapsigargin active time 
disp('Select 2 point on graph to define time point of ER depletion');
ER_dcm = datacursormode(gcf);
set(ER_dcm,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');

disp('Click TWO points on the plot, then press ENTER in the command window...');

pause   % waits for user to finish selecting

% Get cursor info
ER_info = getCursorInfo(ER_dcm);

ER_dep = ER_info(2).DataIndex;
ER_dep(end+1) = ER_info(1).DataIndex;


neuronts_thap = neuronts(:,ER_dep(1):ER_dep(2));

figure(14),clf
plot(timevec(ER_dep(1):ER_dep(2)),neuronts_thap)

auc = zeros(nROIs,1); %preallocate AUC array

for i=1:nROIs
   
    [~,tg_locs] = findpeaks(neuronts_thap(i,:),'MinPeakProminence',1);

    if ~isnan(tg_locs) 
          tg_int = find(neuronts_thap(i,1:tg_locs(1))<=0,1, 'last');
          tg_ov = tg_locs(end) + find(neuronts_thap(i,tg_locs(end):end)<=0,1, 'first');
          auc(i) = trapz(neuronts_thap(i,tg_int:end));
      
    end    

end

auc(auc == 0) = nan; %exclude inactive cells

%visualize as bar - quick check
figure(15),clf
bar(auc)

auc_mean = nanmean(auc);

%% SOCE slope
%define SOCE slope active time 
disp('Select 2 point on graph to define slope');
SOCE_dcm = datacursormode(gcf);
set(SOCE_dcm,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');

disp('Click TWO points on the plot, then press ENTER in the command window...');

pause   % waits for user to finish selecting

% Get cursor info
SOCE_info = getCursorInfo(SOCE_dcm);

SOCE = SOCE_info(2).DataIndex;
SOCE(end+1) = SOCE_info(1).DataIndex;


neuronts_slpe = neuronts(:,SOCE(2))-neuronts(:,SOCE(1));
timevec_slpe = timevec(SOCE(2))-timevec(SOCE(1));
slpe = neuronts_slpe / timevec_slpe;
slpe_mean = mean(slpe);

%% VOCC+SOCE AUC - get start and end time point

% Enable data cursor mode
dcm = datacursormode(gcf);
set(dcm,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');

disp('Click TWO points on the plot, then press ENTER in the command window...');

pause   % waits for user to finish selecting

% Get cursor info
c_info = getCursorInfo(dcm);

% Check number of selected points
if length(c_info) >= 2
    bg_voso  = c_info(2).DataIndex;   % X of first selected point
    end_voso = c_info(1).DataIndex;   % X of second selected point
   

    fprintf('bg_voso  = %d \n', bg_voso);
    fprintf('end_voso = %d \n', end_voso);
else
    warning('Please select at least TWO points.')
    bg_voso  = [];
    end_voso = [];
end

%% VOCC+SOCE AUC - Compute 
neuront_voso = neuronts(:,bg_voso:end_voso)-mean(neuronts(:,bg_voso-20:bg_voso-1),2);

figure(16),clf
plot(timevec(bg_voso:end_voso),neuront_voso)
set(gca,'xlim',timevec([bg_voso end_voso]))

auc_voso = trapz(neuront_voso,2);

figure(17),clf
bar(auc_voso)


%% (Option) Save data
[~, name, ~] = fileparts(filename);
datFileName = ['DATA' name ];
save(datFileName, '-v7.3');

%% end