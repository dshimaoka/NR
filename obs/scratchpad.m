

%% DS OFFICE
% params.paths = '/home/daisuke/Documents/git/NR/Z:/Shared/MarmosetData/2024/07/30/';
% params.files = 'test.nystagmusRivalry.164731.mat';

%% DOWNSTAIRS
params.paths = '/mnt/dshi0006_market/MarmosetData/2024/07/30';
params.files = 'test.nystagmusRivalry.185754.mat'; %grating
patchType = 'grating';

% params.files = 'test.nystagmusRivalry.184310.mat'; %rdp
% patchType = 'rdp';

%data = edfmex( '/mnt/syncitium/Daisuke/temp/marmolab/2024/07/30/test.nystagmusRivalry.185754.edf');

d = marmodata.mdbase('path',params.paths,'file',params.files,'loadArgs',{'loadEye',true});

color1 = d.meta.patch1.color('time',Inf).data;
color2 = d.meta.patch2.color('time',Inf).data;
redFirst = d.meta.patch1.redFirst('time',Inf).data;

%% eye trajectory of each trial
%stimDir = get(c.gabor1.prms.orientation,'atTrialTime',inf);
% stimDir = get(d.meta.patch1.direction, 'atTrialTime', Inf);
stimDir = d.meta.patch1.direction('time',Inf).data;
congruent = d.meta.patch2.congruent('time',Inf).data;
%directionPolarity1 = d.meta.patch1.directionPolarity('time',Inf).data;
%directionPolarity2 = d.meta.patch2.directionPolarity('time',Inf).data;
if strcmp(patchType,'grating')
    radius = unique(d.meta.patch1.sigma('time',Inf).data);
elseif strcmp(patchType,'rdp')
    radius = unique(d.meta.patch1.maxRadius('time',Inf).data);
end
%TODO: 
% use the same name as much as possible betweeen RDP and grating
% stop recording unfactorized parameters
stimOn1 = d.meta.patch1.on('time',Inf).data; %ms
stimOn2 = d.meta.patch2.on('time',Inf).data; %ms

%% retrieve keyPress reporting perceptual switches
% ALL NAN so far
%from https://github.com/nPrice-lab/pursuit/blob/e088cb8747ed8e68bb5c43d69dc8adc14b552e6f/pursuit2D_Merge.m#L45
[time,trial,frame,keyTmp] = d.meta.keypress.keyIx('time',Inf);
key = cell2mat(keyTmp);
ignoreTrial = isnan(key);

keyPressTime = cell(d.numTrials, 1);
for itr = 1:d.numTrials
    keyPressTime{itr} = time(trial == itr);
end

%% successful or failed fixation to the end of the trial
% to be housed in analysis.NR
[~,trial,~,state] = d.meta.fixbhv.state;
outCome = false([d.numTrials,1]);
ix = strcmpi(state,'SUCCESS');
outCome(trial(ix)) = true;
iy = strcmpi(state,'FAIL');
outCome(trial(iy)) = false;

trialEndTime = zeros(d.numTrials,1);
for itr = 1:d.numTrials
    trialEndTime(itr) =  d.eye(itr).t(end); %s
end


%% judge when the trial aborted
numAbort(1) = sum((1e3*trialEndTime <= stimOn1) .* (outCome == false)); %before 1st patch 
numAbort(2) = sum((1e3*trialEndTime > stimOn1) .* (1e3*trialEndTime <= stimOn2) .* (outCome == false)); %between 1st and 2nd patches 
numAbort(3) = sum((1e3*trialEndTime > stimOn2) .* (outCome == false)); %after 2nd patch
sprintf('# aborted trials before 1st patch: %d, between patches: %d, after 2nd patch: %d', numAbort(1), numAbort(2), numAbort(3))


%% eye direction during 1st and 2nd patch


%% remove blinks then saccades then show traces in x&y
for itr = 1:d.numTrials
    tmp = d.eye(itr).rmBlinks;
    eye_rmSaccades(itr) = tmp.rmSaccades;
    close all
end

%% judge if and when the reversal of eye direction occurs
switchTime = cell(d.numTrials, 1);
vel_stimDir = cell(d.numTrials, 1);
t_r = cell(d.numTrials, 1);
for itr = 1:d.numTrials
  
    tidx = find((d.eye(itr).t>0) & (d.eye(itr).t<max(d.eye(itr).t)));

    t = d.eye(itr).t(tidx); %s

    for ii = 1:2
        switch ii
            case 1
                x = eye_rmSaccades(itr).x(tidx);
                signal = x;
            case 2
                y = eye_rmSaccades(itr).y(tidx);
                signal = y;
        end

        okIdx = ~isnan(signal);
        signal = interp1(find(okIdx), signal(okIdx), 1:numel(signal), 'nearest','extrap')'; %extrapolate r to remove NaNs

        order = 3;
        fs = 1/median(diff(t));
        cutoffFreq = 1;%Hz
        Wn = cutoffFreq/(fs/2);
        [b,a]=butter(order, Wn, 'low');

        signal_c = filtfilt(b,a,cat(1,flipud(signal), ...
            signal, flipud(signal)));
        ntotFrames = numel(signal);
        signal_f = signal_c(ntotFrames+1:2*ntotFrames); %filtered distance from fixation cross
        
        fs_r = 100; %resampled frequency [Hz]
        t_r{itr} = t(1):1/fs_r:t(end);
        signal_rf = interp1(t, signal_f,t_r{itr})';
        switch ii
            case 1
                x_f = signal_f; %filtered distance from fixation cross
                x_rf = signal_rf; %resampled
                vx_rf = gradient(x_rf, 1/fs_r);
            case 2
                y_f = signal_f; %filtered distance from fixation cross
                y_rf = signal_rf; %resampled
                vy_rf = gradient(y_rf, 1/fs_r);
        end
    end
    theta = atan2(vy_rf, vx_rf);
    vel_stimDir{itr} = sqrt(vx_rf.^2+vy_rf.^2) .* cos(theta - pi/180*stimDir(itr));

[pks, peakTidx_c] =     findpeaks(-abs(vel_stimDir{itr}), 'MinPeakProminence',3);
peakTidx = [];
dt_ba = 0.1;%[s] %time points to judge polarity change before and after the peak
for jj = 1:numel(pks)
    if pks(jj) < 0.1 && peakTidx_c(jj) > dt_ba*fs_r+1 && peakTidx_c(jj) < numel(t_r{itr})-dt_ba*fs_r && ...
            vel_stimDir{itr}(peakTidx_c(jj)-dt_ba*fs_r)*vel_stimDir{itr}(peakTidx_c(jj)+dt_ba*fs_r)<0
        peakTidx = [peakTidx peakTidx_c(jj)];
    end
end
    switchTime{itr} = 1e3*t_r{itr}(peakTidx);
end
switched = cellfun(@(x)~isempty(x), switchTime);

%% difference in direction between stimulus and bhv


%% gain of OKN



%% show bhv performance by stimulus direction
stimDirList = unique(stimDir);
nCompleteTrials = zeros(numel(stimDirList), 2, 2);
nSwitchedTrials = zeros(numel(stimDirList), 2, 2);
for idir = 1:numel(stimDirList)
    for ifirstRed = 0:1
        for icongruent = 0:1
            completeTrials = (stimDir == stimDirList(idir)) .* ...
                (redFirst==ifirstRed) .* (outCome==true) .*  (congruent == icongruent);
            nCompleteTrials(idir, ifirstRed+1, icongruent+1) = sum( completeTrials); %supposed to be all d.numTrials/d.numConds

            nSwitchedTrials(idir, ifirstRed+1, icongruent+1) = sum((switched==1) .* completeTrials);
        end
    end
end
figure;
plot(stimDirList, nSwitchedTrials(:,1,1)./nCompleteTrials(:,1,1), '-s','Color','b','MarkerFaceColor','b'); hold on
plot(stimDirList, nSwitchedTrials(:,1,2)./nCompleteTrials(:,1,2), '-o','Color','b');
plot(stimDirList, nSwitchedTrials(:,2,1)./nCompleteTrials(:,2,1), '-s','Color','r','MarkerFaceColor','r'); 
plot(stimDirList, nSwitchedTrials(:,2,2)./nCompleteTrials(:,2,2), '-o','Color','r');
axis padded
legend('blueFirst incongruent', 'blueFirst congruent','redFirst incongruent', 'redFirst congruent','Location','eastoutside');
xlabel('stimulus direction [deg]'); ylabel('eye switch rate');
screen2png(['nSwitchedTrials_' params.files '.png']);


%% sannity check of a given trial
winSize = 10;
for itr = 1:d.numTrials

    tidx = find((d.eye(itr).t>0) & (d.eye(itr).t<max(d.eye(itr).t)));
    t = d.eye(itr).t(tidx); %s

    [~, stimOn1tidx] = min(abs(1e3*t - stimOn1(itr)));
    [~, stimOn2tidx] = min(abs(1e3*t - stimOn2(itr)));

    for ii = 1:2
        switch ii
            case 1
                x = d.eye(itr).x(tidx);
                x_rm = eye_rmSaccades(itr).x(tidx);
            case 2
                y = d.eye(itr).y(tidx);
                y_rm = eye_rmSaccades(itr).y(tidx);
        end
    end
    r = sqrt(x_rm.^2+y_rm.^2); %distance from the fixation point
    okIdx = ~isnan(r);
    r = interp1(find(okIdx), r(okIdx), 1:numel(r), 'nearest','extrap')'; %extrapolate r to remove NaNs

    order = 3;
    fs = 1/median(diff(t));
    cutoffFreq = 1;%Hz
    Wn = cutoffFreq/(fs/2);
    [b,a]=butter(order, Wn, 'low');
    signal_c = filtfilt(b,a,cat(1,flipud(r), ...
        r, flipud(r)));
    ntotFrames = numel(r);
    r_f = signal_c(ntotFrames+1:2*ntotFrames); %filtered distance from fixation cross

    figure('position', [1 1 750 450]);
    subplot(221);
    plot(1e3*t, x, 1e3*t, y); legend('x','y');
    if congruent(itr) == 1
        thisColor = 'k';
    else
        thisColor = 'r';
    end
    vline(stimOn2(itr),gca,'-',thisColor);
    vbox(stimOn1(itr),1e3*t(end),gca,[redFirst(itr) 0.0 1-redFirst(itr) .1])
    vbox(stimOn2(itr),1e3*t(end),gca,[1-redFirst(itr) 0.0 redFirst(itr) .1])
    vline(switchTime{itr}, gca, '--','g');
    title(sprintf('tr:%d', itr))
    xlabel('time [ms]'); ylabel('raw eye position [deg]')
    if outCome(itr)==0; set(gca,'xcolor','r','ycolor','r');end

    subplot(222);
    plot(1e3*t, x_rm, 1e3*t, y_rm); 
    vline(stimOn2(itr),gca,'-',thisColor);
    vbox(stimOn1(itr),1e3*t(end),gca,[redFirst(itr) 0.0 1-redFirst(itr) .1])
    vbox(stimOn2(itr),1e3*t(end),gca,[1-redFirst(itr) 0.0 redFirst(itr) .1])
    vline(switchTime{itr}, gca, '--','g'); 
    xlabel('time [ms]'); ylabel('after removal of saccade [deg]');
 if outCome(itr)==0; set(gca,'xcolor','r','ycolor','r');end

    subplot(223);
    %plot(1e3*t, r, 'k',1e3*t, r_f, 'g');
    plot(1e3*t_r{itr}, vel_stimDir{itr});
    vline(stimOn2(itr),gca,'-',thisColor);
    vbox(stimOn1(itr),1e3*t(end),gca,[redFirst(itr) 0.0 1-redFirst(itr) .1])
    vbox(stimOn2(itr),1e3*t(end),gca,[1-redFirst(itr) 0.0 redFirst(itr) .1])
    hline(0);
    vline(switchTime{itr}, gca, '--','g');
    xlabel('time [ms]');
    ylabel('velocity along 1st stim');
 if outCome(itr)==0; set(gca,'xcolor','r','ycolor','r');end

    subplot(224)
    quiver(0,0, winSize*cos(stimDir(itr)*pi/180), winSize*sin(stimDir(itr)*pi/180),'linewidth',1); hold on;
    try
        scatter(x_rm(stimOn1tidx:end) - x_rm(stimOn1tidx), y_rm(stimOn1tidx:end) - y_rm(stimOn1tidx), 1,t(stimOn1tidx:end));
        plot(x_rm(stimOn2tidx)- x_rm(stimOn1tidx), y_rm(stimOn2tidx)- y_rm(stimOn1tidx), 'o', 'Color',thisColor);
    end
    xlim([-winSize winSize]);  ylim([-winSize winSize]);
    vline(0,gca,'--',[.5 .5 .5 .5]); hline(0,gca,'--',[.5 .5 .5 .5]);
    xlabel('x [deg]'); ylabel('y [deg]');
    axis square;
    if outCome(itr)==0; set(gca,'xcolor','r','ycolor','r');end

    screen2png(['eyeSummary_' params.files '_tr' num2str(itr) '.png']);
    close
end
