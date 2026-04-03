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