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


%%
%select ER deplete start point
star_ER_dcm = datacursormode(gcf);
set(star_ER_dcm,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');
pause   % waits for user to finish selecting
star_ER_info = getCursorInfo(star_ER_dcm);

star_ER_dep = star_ER_info(1).DataIndex;

meanNorm_ER = meanNorm(star_ER_dep-20:star_ER_dep+240);
semNorm_ER = semNorm(star_ER_dep-20:star_ER_dep+240);
timevecplot_ER = (1:1:size(meanNorm_ER))/2; 

%%
figure(1000)
fill([timevecplot_ER'; flipud(timevecplot_ER')], ...
     [meanNorm_ER+semNorm_ER; flipud(meanNorm_ER-semNorm_ER)], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(timevecplot_ER, meanNorm_ER, 'b', 'LineWidth', 1.5)
ylabel('\DeltaF/F_0')
xlabel('Time (sec)')
title('\DeltaF/F_0 Trace with SEM')
set(gca, 'xlim', timevecplot_ER([1 80]) , ...
    'ylim', [-1.5 15], ...
    'YTick',(0:2:15), ...
    'fontsize',12)
set(gca,'plotboxaspectratio',[1 4 1])
hold on