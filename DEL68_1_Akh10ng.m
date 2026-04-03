%% Visualize time series data from table
% AI generated Load Excel data
% Visualization by Mike X Cohen, sincxpress.com (source code)
% Slightly modification in normalization step by Thanh
close all, clear, clc
%% Load Excel time series data with multiple ROIs
filename = ['DEL68-Rsc68-3-Akh1ng.csv'];

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
evnt= {'Akh 1 ng'};


%% Visualization

% Plot an example trace
figure(4), clf
%rprest = randi(nROIs,1);
j= 5 ;
plot(timevec, neuronts(j,:))
ylabel('Brightness (a.u.)')
set(gca, 'xlim', timevec([1 end]), ...
    'ylim', [0 1000])
title(sprintf('Fluorescence of ROIs %d (normalized)', j))
box off

% Show heatmap
figure(5), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),neuronts,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5])
set(gca, 'xlim', timevec([1 end]))
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
set(gca, 'xlim', timevec([1 end]), ...
    'YLim',[0 1500])

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
set(gca, 'xlim', timevec([1 end]))
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
set(gca, 'xlim', timevec([1 end]))
%% Normalized to tomato
neuronts = bsxfun(@rdivide,neuronts,mean(refneuronts(:,40:60),2));
%neuronts = bsxfun(@rdivide,neuronts,refneuronts);

base_ca = mean(neuronts(:,24:marks(1)),2);
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


%% Convert to dF/F

neuronts = bsxfun(@rdivide, neuronts, mean(neuronts(:,marks(1)-20:marks(1)-1), 2))-1;

%% dF/F visualization
figure(11), clf
surf(timevec,linspace(0.5,nROIs+0.5,nROIs),neuronts,'LineStyle','none')
set(gca,'ylim',[0.5 nROIs+0.5], ...
    'ytick',[1 nROIs],...
        'xlim', timevec([1 end]), ...
    'clim',[-1 10], ...
    'YDir','reverse', ...
    'fontsize',20, ...
    'fontname','Arial Narrow')
xline(timevec(marks),'LineWidth',1.5,'Color',[1 1 1],'LineStyle','--')
text(timevec(marks(1)),-0.04,evnt,"HorizontalAlignment","left",'FontSize',20)
colorbar
axis square
box off
grid off
view(2)

figure(12), clf
set(gcf,'Units', 'inches')
set(gcf,'Position', [1 1 15 5])
plot(timevec, neuronts)
ylabel('\DeltaF/F_0')
xlabel('Time (min.)')
title('Fluorescence of all ROIs (normalized)', ...
    'FontWeight','normal', ...
    'fontsize',22, ...
    'FontName','Arial Narrow')
set(gca,'xlim', timevec([1 end]), ...
    'ylim', [-1.5 12], ...
    'YTick',(0:2:10), ...
    'TickDir','none', ...
    'Fontsize',20, ...
    'Fontname','Arial Narrow')
box off
%drawing 1 lines to indicate drugs treatment on top of graph
line([timevec(marks(1)) timevec(end)],[10.5 10.5],'color','k','Linewidth',2)
text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),11.5,evnt, ...
    'HorizontalAlignment','center','FontSize',20)

%% (Optional) all cell time tracing
for i=1:nROIs
    figure(100+i), clf
    set(gcf,'Units', 'inches')
    set(gcf,'Position', [1 1 15 5])
    plot(timevec, neuronts(i,:),'-k','LineWidth',1)
    ylabel('\DeltaF/F_0','FontSize',18)
    xlabel('Time (min.)','FontSize',18)
    title(sprintf('Fluorescence of ROIs %d (normalized)', i))
    set(gca, 'xlim', timevec([1 end]),'FontSize',16)
    set(gca, 'ylim', [-0.8 10],'FontSize',16)
    xline(timevec([marks]),'--b',evnt,'fontsize',13,'FontWeight','bold','LineWidth',1)

end

%% (Option) 1 cell tracing
    j =5; %define cell
    figure(200+j), clf
    set(gcf,'Units', 'inches')
    set(gcf,'Position', [1 1 15 5])
    plot(timevec, neuronts(j,:),'-k','LineWidth',1)
    ylabel('\DeltaF/F_0','FontSize',18)
    xlabel('Time (min.)','FontSize',18)
    title(sprintf('Fluorescence of ROIs %d (normalized)', j))
    set(gca,'xlim', timevec([1 end]), ...
    'ylim', [-1.5 12], ...
    'YTick',(0:2:10), ...
    'TickDir','none', ...
    'Fontsize',20, ...
    'Fontname','Arial Narrow')
    box off
    %drawing 1 lines to indicate drugs treatment on top of graph
    line([timevec(marks(1)) timevec(end)],[10.5 10.5],'color','k','Linewidth',2)
    text(timevec(marks(1))+0.5*(timevec(end)- timevec(marks(1))),11.5,evnt, ...
        'HorizontalAlignment','center','FontSize',20)
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
subplot(2,1,1)
hold on
fill([timevec; flipud(timevec)], ...
     [meanRaw+semRaw; flipud(meanRaw-semRaw)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(timevec, meanRaw, 'b', 'LineWidth', 1.5)
ylabel('Raw Brightness (a.u.)')
title('Raw Trace with SEM')
set(gca, 'xlim', timevec([1 end]))
xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold')
box off

% Normalized
subplot(2,1,2)
hold on
fill([timevec; flipud(timevec)], ...
     [meanNorm+semNorm; flipud(meanNorm-semNorm)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(timevec, meanNorm, 'r', 'LineWidth', 1.5)
ylabel('\DeltaF/F_0')
xlabel('Time (min)')
title('\DeltaF/F_0 Trace with SEM')
set(gca, 'xlim', timevec([1 end])) %, 'ylim', [-0.2 4])
xline(timevec([marks]),'--b',evnt,'fontsize',11,'FontWeight','bold')
box off


%% peak analyze 100ng
%%
min_height = 1;
min_pros = 1;
min_wds = 0.1;
%dif_thr = 0.001;
Akh_conc = [10];
min_dst = 1;

% Initialize cell arrays
hdr_pks = cell(size(neuronts,1),1);
hdr_locs = cell(size(neuronts,1),1);

for i= 1:size(neuronts,1)
     % Extract signal and time window
    signal_100 = neuronts(i, marks(1):end); %define time point
    time_100 = timevec(marks(1):end);
    %exclude any cell with failed washing
    if mean(neuronts(i,40:marks(1)),2)<1
         % Find peaks with thresholds
        [pks_100,locs_100] = findpeaks(signal_100,time_100, 'MinPeakHeight', min_height, ...
            'MinPeakProminence', min_pros,'MinPeakWidth',min_wds,'MinPeakDistance',min_dst );
         % Keep only the first peak if it's at least 2× taller than others
        %if ~isempty(pks_100)
         %   if pks_100(1) > 1.8 * max(pks_100(2:end))
          %      pks_100 = pks_100(1);
           %     locs_100 = locs_100(1);
            %end
        %end
    %visualize (doublecheck)
    figure(i), clf
    %raw
    subplot(2,1,1)
    plot(time_100,neuronts_raw(i,marks(1):end))
     set(gca, 'xlim',timevec([marks(1) end]),'FontSize',16)

    %df
    subplot(2,1,2)
    plot(time_100,signal_100, 'LineWidth',1)
    hold on
    plot(locs_100,pks_100+0.35, 'rv','MarkerFaceColor','r','Marker','v', ...
               'LineStyle','none')
    set(gca, 'xlim',timevec([marks(1) end]),'FontSize',16)
    set(gca,'ylim',[-1 15],'FontSize',16)
    %Store entire peak array in one cell
    hdr_pks{i} = pks_100;
    hdr_locs{i} = locs_100;
    else
    %do nothing
    end
end

%% one concetration (100 ng Akh)
peak_counts = cellfun(@numel, hdr_pks );
peak_counts(peak_counts == 0) = NaN; 
mean_oscil = nanmean(peak_counts,1 );
peak_ave = cellfun(@mean, hdr_pks);
ave_amp = nanmean(peak_ave,1 );
act_cell = sum(peak_counts>0);
%calculate percentage of activated cell
per_act_cell = act_cell ./size(neuronts,1);
per_act_cell = per_act_cell .* 100;
sem_amp  = std(peak_ave,0,1,'omitmissing') / sqrt(act_cell);
sem_oscil  = std(peak_counts,0,1,'omitmissing') / sqrt(act_cell);

%% get the 1st peaks
%get 1st peaks from previous variable
st_hdr_locs = zeros(nROIs,1);
st_hdr_pks = zeros(nROIs,1);
st_int = zeros(nROIs,1);
st_loss = zeros(nROIs,1);
for k =1:nROIs
    if ~isnan(hdr_locs{k}) %&(hdr_locs{k}(1)<10)
        
        st_hdr_locs(k) = hdr_locs{k}(1);
        st_hdr_pks(k) = hdr_pks{k}(1);
        
        %define int and loss
        %frame > time = 5*12+ x*30 bsc 5 min baseline has 60 frame and then
        %record signal in every 2 secs
        st_int_bline =  prctile(neuronts(k,marks(1):find(timevec==st_hdr_locs(k))),15);
        st_int(k) = timevec(marks(1)+find(neuronts(k,marks(1):find(timevec==st_hdr_locs(k)))<= ...
            st_int_bline,1,'last')); %minutes
        
       if isempty(find(neuronts(k,find(timevec==st_hdr_locs(k)):end)<=0.3))==1;
           st_loss(k) = 35;

       else
           if length(hdr_pks{k}) ==1
              st_loss(k) = timevec(find(timevec==st_hdr_locs(k)) + find(neuronts(k,find(timevec==st_hdr_locs(k)):end) ...
                  <= ...
                  prctile(neuronts(k,find(timevec==st_hdr_locs(k)):end),15),1,'first'));
              st_to_bline = prctile(neuronts(k,find(timevec==st_hdr_locs(k)):end),15);
           else
               st_nd_dy = neuronts(k,find(timevec==st_hdr_locs(k)):find(timevec==hdr_locs{k}(2)));
               st_nd_bline=prctile(st_nd_dy,18);
               st_to_bline = find(st_nd_dy <= st_nd_bline,1,'first');
               st_loss(k) = timevec(find(timevec==st_hdr_locs(k))+st_to_bline);
               
           
           end

           
           figure(400+k),clf
           hold on
           plot(timevec, neuronts(k,:),'-k','LineWidth',1)
           plot(hdr_locs{k},hdr_pks{k}+0.3,'MarkerFaceColor','r','Marker','v', ...
               'LineStyle','none')
           ylabel('\DeltaF/F_0','FontSize',18)
           xlabel('Time (min.)','FontSize',18)
           title(sprintf('Fluorescence of ROIs %d (normalized)', k))
           set(gca, 'xlim', timevec([1 end]),'FontSize',16)
           set(gca, 'ylim', [-0.8 10],'FontSize',16)
           xline(timevec(marks(1)),'--b',evnt,'fontsize',11,'FontWeight','bold')
           rectangle('Position',[st_int(k) st_int_bline st_hdr_locs(k)-st_int(k) ...
               st_hdr_pks(k)+0.35-(st_int_bline)],'FaceColor','y','FaceAlpha',.5,'EdgeColor','none')
           rectangle('Position',[st_hdr_locs(k) st_nd_bline st_loss(k)-st_hdr_locs(k) ...
               st_hdr_pks(k)+0.35-(st_nd_bline)],'FaceColor','m','FaceAlpha',.5,'EdgeColor','none')
           %xline(st_int(k),'g-')
           %yline(st_nd_bline)
           

       end
    
    end
end

st_hdr_pks(st_hdr_pks == 0) = NaN;
st_hdr_locs (st_hdr_locs ==0) = NaN;
%compute 1st peaks high
st_amp = nanmean(st_hdr_pks);
st_loc = nanmean(st_hdr_locs);

%compute raise and decay and duration
st_int(st_int ==0) = NaN;
st_loss(st_loss ==0) = NaN;

st_raise = st_hdr_locs -  st_int;
st_decay = st_loss - st_hdr_locs;
st_dur = st_loss - st_int;

disp('completed 1st peaks');


%% from 2nd peaks to the end
%get peaks amp and locs
rest_locs =cell(nROIs,1);
rest_pks = cell(nROIs, 1);
for l=1:nROIs
    if length(hdr_locs{l}) > 1
        rest_locs{l} = hdr_locs{l}(2:end);
        rest_pks{l} = hdr_pks{l}(2:end);
    end
end

rest_cnt_osc = cellfun(@numel, rest_pks);
rest_cnt_osc(rest_cnt_osc ==0) = NaN;
rest_cnt_amp = cellfun(@mean, rest_pks);

%% (Option) Save data
[~, name, ~] = fileparts(filename);
datFileName = ['DATA' name ];
save(datFileName, '-v7.3');

%% end