function d = summary(file)

 set(0,'DefaultFigureVisible','off');
% if nargin < 2
%     saveServer = '/mnt/dshi0006_market/MarmosetAnalysis';
% end

[saveDir, thisFile] = fileparts(file);

% d = marmodata.mdbase('path',params.paths,'file',params.files,'loadArgs',{'loadEye',true});
d = analysis.NR('file',file,'loadArgs',{'loadEye',true}); %loads eye data into 'd' 
disp('data loaded!') 

disp([d.patchType ': ' d.file{:}]);

%% aborted trials
d.getNumAbort;

%% complete but invalid due to eye movements
d.getInvalidBhvTrials;

%% eye switch latency FS v PA
d.statEyeSwitch;

%% eye - key response consistency
d.eyeKeyConsistency;

%% eye trace per trial
disp('Plotting single trial eye traces');
for itr = 1:d.numTrials
    d.plotSingleTrialEye(itr);
    screen2png(fullfile(saveDir, ['eyeSummary_' thisFile '_tr' num2str(itr) '.png']));
    close;
end
disp('done');

%% eye switch by stimulus direction
d.plotSwitchByPatchDir;
screen2png(fullfile(saveDir,['nSwitchedTrials_' thisFile '.png']));
close

 set(0,'DefaultFigureVisible','on');