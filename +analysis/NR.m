classdef NR < marmodata.mdbase % vgsaccade.vgsaccade
    % base class for analysis of the ocular following paradigm with motion
    % cloud stim

    properties
        % parameters for eye movement switch
        cutoffFreq = 1;%low pass filtering freq for detection of eye movement switch [Hz]
        dt_ba = 0.1;%[s] %time points to judge polarity change before and after the peak
        fs_r = 100; %resampled frequency [Hz]
        prominenceThresh = 2; %threshold for velocity in detecting eye movement switch[deg/s]
        accThresh = 2500; %5000; %threshold for saccade detection in findSaccades
    end

    properties
        % stimulus property
        redFirst; %whether the red patch is presented first (1) or second (0)
        stim1Dir; %direction of a first patch [deg]
        congruent; %whether the two patches are same directions
        patchType; %kinds of patches ('grating' or 'rdp')
        radius; %radius of patches in [deg]
        stimDirList; %list of stimulus directions
        stimSpeed; %stimulus speed [deg/s]

        % stimulus timing
        stim1On; %onset time of 1st patch onset [ms]
        stim2On; %onset time of 2nd patch onset [ms]

        % eye data
        eye_rm; %eye data after removing blinks and saccades

        % behavioural timing
        switchTime; %time when eye movement switches [ms]

        % screening trials
        complete; %whether a subject continued fixation till the end of a trial


    end

    methods  (Access = public)
        function d = NR(varargin)
            fprintf(1,'@NR.NR(): nargin = %i.\n', nargin);
            d@marmodata.mdbase(varargin{:}); % base class constructor
            d.redFirst = getRedFirst(d);
            d.complete = getComplete(d);
            d.stim1Dir = getStim1Dir(d);
            d.stimDirList = getStimDirList(d);
            d.congruent = getCongruent(d);
            d.patchType = getPatchType(d);
            d.radius = getRadius(d);
            d.stim1On = getStim1On(d);
            d.stim2On = getStim2On(d);
            d.stimSpeed = getStimSpeed(d);
            d.eye_rm = rmBlinkSaccade(d);
            d.switchTime = getSwitchTime(d);
        end

        % function tDur = getTDur(d)
        % end
        % function nRepPerCond = getNRepPerCond(d)
        % end
        % function rewardVol = getRewardVol(d)
        % end
        % function conditionSwitch = getConditionSwitch(d)
        % end
        % function SOARange = getSOARange(d)
        % end

        function speed = getStimSpeed(d)
            speed = unique(d.meta.patch1.speed('time',Inf).data);
        end

        function v = getComplete(d)
            [~,trial,~,state] = d.meta.fixbhv.state;
            v = false([d.numTrials,1]);
            ix = strcmpi(state,'SUCCESS');
            v(trial(ix)) = true;
            iy = strcmpi(state,'FAIL');
            v(trial(iy)) = false;
        end

        function redFirst = getRedFirst(d)
            redFirst = d.meta.patch1.redFirst('time',Inf).data;
        end

        function stim1Dir = getStim1Dir(d)
            stim1Dir = d.meta.patch1.direction('time',Inf).data;
        end

        function stimDirList = getStimDirList(d)
            stimDirList = unique(d.stim1Dir);
        end

        function congruent = getCongruent(d)
            congruent = d.meta.patch2.congruent('time',Inf).data;
        end

        function radius = getRadius(d)
            if strcmp(d.patchType,'grating')
                radius = unique(d.meta.patch1.sigma('time',Inf).data);
            elseif strcmp(d.patchType,'rdp')
                radius = unique(d.meta.patch1.maxRadius('time',Inf).data);
            end
        end

        function patchType = getPatchType(d)
            if sum(strcmp(d.meta.patch1.show, 'patch1.sigma'))
                patchType = 'grating';
            elseif sum(strcmp(d.meta.patch1.show, 'patch1.maxRadius'))
                patchType = 'rdp';
            else
                patchType = 'unknown';
            end
        end

        function stim1On = getStim1On(d)
            stim1On = d.meta.patch1.on('time',Inf).data; %ms
        end

        function stim2On = getStim2On(d)
            stim2On = d.meta.patch2.on('time',Inf).data; %ms
        end

        function keyPressTime = getKeyPressTime(d)
            %get time of key press(es)
            % TOBE CONFIRMED
            [time,trial,frame,keyTmp] = d.meta.keypress.keyIx('time',Inf);
            key = cell2mat(keyTmp);
            ignoreTrial = isnan(key);

            keyPressTime = cell(d.numTrials, 1);
            for itr = 1:d.numTrials
                keyPressTime{itr} = time(trial == itr);
            end
        end

        function trialEndTime = getTrialEndTime(d)
            trialEndTime = zeros(d.numTrials,1);
            for itr = 1:d.numTrials
                trialEndTime(itr) =  d.eye(itr).t(end); %s
            end
        end

        function numAbort = getNumAbort(d)
            trialEndTime = d.getTrialEndTime;
            numAbort(1) = sum((1e3*trialEndTime <= d.stim1On) .* (d.complete == false)); %before 1st patch
            numAbort(2) = sum((1e3*trialEndTime > d.stim1On) .* (1e3*trialEndTime <= d.stim2On) .* (d.complete == false)); %between 1st and 2nd patches
            numAbort(3) = sum((1e3*trialEndTime > d.stim2On) .* (d.complete == false)); %after 2nd patch
            sprintf('# aborted trials before 1st patch: %d, between patches: %d, after 2nd patch: %d', numAbort(1), numAbort(2), numAbort(3))
        end

        function eye_rm = rmBlinkSaccade(d)
            for itr = 1:d.numTrials
                tmp = d.eye(itr).rmBlinks('dt',median(diff(d.eye(itr).t)), 'duration', 0.01, 'debug', false); %marmodata/+marmodata/@eye/rmBlinks.m
                eye_rm(itr,1) = tmp.rmSaccades('debug',false,'sargs',{'accthresh', d.accThresh}); %marmodata/+marmodata/@eye/rmSaccades.m
                %cf. fitKernel/selectSaccades
                %close all
            end
        end

        function [switchTime, vel_stimDir, avgGain_stimDir] = getSwitchTime(d)
            switchTime = cell(d.numTrials, 1); %ms
            vel_stimDir = cell(d.numTrials, 1); % deg/s
            avgGain_stimDir = nan(d.numTrials, 1);
            for itr = 1:d.numTrials

                tidx = find((d.eye(itr).t>0) & (d.eye(itr).t<max(d.eye(itr).t)));

                t = 1e3*d.eye(itr).t(tidx); %s

                for ii = 1:2
                    switch ii
                        case 1
                            x = d.eye_rm(itr).x(tidx);
                            signal = x;
                        case 2
                            y = d.eye_rm(itr).y(tidx);
                            signal = y;
                    end

                    okIdx = ~isnan(signal);
                    signal = interp1(find(okIdx), signal(okIdx), 1:numel(signal), 'nearest','extrap')'; %extrapolate r to remove NaNs

                    signal_f = analysis.src.lowpassFilter(signal, 1e-3*t, d.cutoffFreq);

                    t_r = t(1):1e3*1/d.fs_r:t(end);
                    signal_rf = interp1(t, signal_f,t_r)';
                    switch ii
                        case 1
                            x_f = signal_f; %filtered distance from fixation cross
                            x_rf = signal_rf; %resampled
                            vx_rf = gradient(x_rf, 1/d.fs_r);
                        case 2
                            y_f = signal_f; %filtered distance from fixation cross
                            y_rf = signal_rf; %resampled
                            vy_rf = gradient(y_rf, 1/d.fs_r);
                    end
                end
                theta = atan2(vy_rf, vx_rf);
                vel_stimDir{itr} = sqrt(vx_rf.^2+vy_rf.^2) .* cos(theta - pi/180*d.stim1Dir(itr));

                [pks, peakTidx_c] =     findpeaks(-abs(vel_stimDir{itr}), 'MinPeakProminence',d.prominenceThresh);
                peakTidx = [];
                for jj = 1:numel(pks)
                    if pks(jj) < 0.1 && peakTidx_c(jj) > d.dt_ba*d.fs_r+1 && peakTidx_c(jj) < numel(t_r)-d.dt_ba*d.fs_r && ...
                            vel_stimDir{itr}(peakTidx_c(jj)-d.dt_ba*d.fs_r)*vel_stimDir{itr}(peakTidx_c(jj)+d.dt_ba*d.fs_r)<0
                        peakTidx = [peakTidx peakTidx_c(jj)];
                    end
                end
                switchTime{itr} = t_r(peakTidx);
                avgGain_stimDir(itr) = median(abs(vel_stimDir{itr}(t_r>d.stim1On(itr))))/d.stimSpeed;
            end
            
        end

        %% functions for visualization
        function fig = plotSwitchByStimDir(d)

            switched = cellfun(@(x)~isempty(x), d.switchTime);

            nCompleteTrials = zeros(numel(d.stimDirList), 2, 2);
            nSwitchedTrials = zeros(numel(d.stimDirList), 2, 2);
            for idir = 1:numel(d.stimDirList)
                for ifirstRed = 0:1
                    for icongruent = 0:1
                        completeTrials = (d.stim1Dir == d.stimDirList(idir)) .* ...
                            (d.redFirst==ifirstRed) .* (d.complete==true) .*  (d.congruent == icongruent);
                        nCompleteTrials(idir, ifirstRed+1, icongruent+1) = sum(completeTrials); %supposed to be all d.numTrials/d.numConds

                        nSwitchedTrials(idir, ifirstRed+1, icongruent+1) = sum((switched==1) .* completeTrials);
                    end
                end
            end
            fig = figure;
            plot(d.stimDirList, nSwitchedTrials(:,1,1)./nCompleteTrials(:,1,1), '-s','Color','b','MarkerFaceColor','b'); hold on
            plot(d.stimDirList, nSwitchedTrials(:,1,2)./nCompleteTrials(:,1,2), '-o','Color','b');
            plot(d.stimDirList, nSwitchedTrials(:,2,1)./nCompleteTrials(:,2,1), '-s','Color','r','MarkerFaceColor','r');
            plot(d.stimDirList, nSwitchedTrials(:,2,2)./nCompleteTrials(:,2,2), '-o','Color','r');
            axis padded
            legend('blueFirst incongruent', 'blueFirst congruent','redFirst incongruent', 'redFirst congruent','Location','eastoutside');
            xlabel('stimulus direction [deg]'); ylabel('eye switch rate');
        end

        function fig = plotSingleTrialEye(d,itr)
            winSize = 10;

            [switchTime, vel_stimDir] = d.getSwitchTime;

            tidx = find((d.eye(itr).t>0) & (d.eye(itr).t<max(d.eye(itr).t)));
            t = 1e3*d.eye(itr).t(tidx); %s
            t_r = t(1):1e3/d.fs_r:t(end);

            [~, stim1Ontidx] = min(abs(t - d.stim1On(itr)));
            [~, stim2Ontidx] = min(abs(t - d.stim2On(itr)));

            for ii = 1:2
                switch ii
                    case 1
                        x = d.eye(itr).x(tidx);
                        x_rm = d.eye_rm(itr).x(tidx);
                    case 2
                        y = d.eye(itr).y(tidx);
                        y_rm = d.eye_rm(itr).y(tidx);
                end
            end

            fig = figure('position', [1 1 750 450]);
            subplot(221);
            plot(t, x,t, y); legend('x','y');
            if d.congruent(itr) == 1
                thisColor = 'k';
            else
                thisColor = 'r';
            end
            vline(d.stim2On(itr),gca,'-',thisColor);
            vbox(d.stim1On(itr),t(end),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.stim2On(itr),t(end),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            vline(switchTime{itr}, gca, '--','g');
            title(sprintf('tr:%d', itr))
            xlabel('time [ms]'); ylabel('raw eye position [deg]')
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(222);
            plot(t, x_rm, t, y_rm);
            vline(d.stim2On(itr),gca,'-',thisColor);
            vbox(d.stim1On(itr),t(end),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.stim2On(itr),t(end),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            vline(switchTime{itr}, gca, '--','g');
            xlabel('time [ms]'); ylabel('after removal of saccade [deg]');
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(223);
            %plot(1e3*t, r, 'k',1e3*t, r_f, 'g');
            plot(t_r, vel_stimDir{itr});
            vline(d.stim2On(itr),gca,'-',thisColor);
            vbox(d.stim1On(itr),t(end),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.stim2On(itr),t(end),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            hline(0);
            vline(switchTime{itr}, gca, '--','g');
            xlabel('time [ms]');
            ylabel('velocity along 1st stim');
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(224)
            quiver(0,0, winSize*cos(d.stim1Dir(itr)*pi/180), winSize*sin(d.stim1Dir(itr)*pi/180),'linewidth',1); hold on;
            try
                scatter(x_rm(stim1Ontidx:end) - x_rm(stim1Ontidx), y_rm(stim1Ontidx:end) - y_rm(stim1Ontidx), 1,t(stim1Ontidx:end));
                plot(x_rm(stim2Ontidx)- x_rm(stim1Ontidx), y_rm(stim2Ontidx)- y_rm(stim1Ontidx), 'o', 'Color',thisColor);
            end
            colormap("cool");
            xlim([-winSize winSize]);  ylim([-winSize winSize]);
            vline(0,gca,'--',[.5 .5 .5 .5]); hline(0,gca,'--',[.5 .5 .5 .5]);
            xlabel('x [deg]'); ylabel('y [deg]');
            axis square;
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end
        end
    end
end