classdef NR < marmodata.mdbase % vgsaccade.vgsaccade
    % base class for analysis of the ocular following paradigm with motion
    % cloud stim

    properties
        %parameters for saccade detection
        accThresh = 2500; %5000; %threshold for saccade detection in findSaccades

        % parameters for eye movement switch
        cutoffFreq = 2;%1;%low pass filtering freq for detection of eye movement switch [Hz]
        dt_ba = 0.1;%[s] %time points to judge polarity change before and after the peak
        fs_r = 100; %resampled frequency [Hz]
        %prominenceThresh = 2; %threshold for velocity in detecting eye movement switch[deg/s]

        % detection of invalid behavioural trials
        gainThresh = 0.05;
        angdiffThresh = 40; %[deg]
    end

    properties

        fixDuration; %duration of initial fixation [ms]

        % stimulus property
        redFirst; %whether the red patch is presented first (1) or second (0)
        patch1Dir; %direction of first patch [deg]
        patch2Dir; %direction of 2nd patch [deg]
        congruent; %whether the two patches are same directions
        patchType; %kinds of patches ('grating' or 'rdp')
        radius; %radius of patches in [deg]
        patchDirList; %list of stimulus directions
        patchSpeed; %stimulus speed [deg/s]
        tDur; %duration of patch1+patch2 [ms]
        conditionSwitchList; %list of stimulus conditions 0: FS, 1: PA, 2: congruent
        conditionSwitch; %condition swith of each trial
       SOA; %time difference between patches 1 and 2 [ms]

        rewardVol;

        % stimulus timing of each trial
        patch1Start; %onset time of 1st patch startTime [ms]
        patch2Start; %onset time of 2nd patch startTime [ms]
        patch1Stop; %onset time of 1st patch stopTime [ms]
        patch2Stop; %onset time of 2nd patch stopTime [ms]

        % eye data
        eye_rm; %eye data after removing blinks and saccades

        % behavioural timing
        switchTime; %time when eye movement switches [ms]
        keyPressTime; %time when a key is perssed first time in a trial [ms]
        % screening trials
        complete; %whether a subject continued fixation till the end of a trial


    end

    methods  (Access = public)
        function d = NR(varargin)
            fprintf(1,'@NR.NR(): nargin = %i.\n', nargin);
            d@marmodata.mdbase(varargin{:}); % base class constructor
            d.redFirst = getRedFirst(d);
            d.congruent = getCongruent(d);
            d.patchType = getPatchType(d);
            d.complete = getComplete(d);
            d.patch1Dir = getPatch1Dir(d);
            d.patch2Dir = getPatch2Dir(d);
            d.patchDirList = getPatchDirList(d);
            d.radius = getRadius(d);
            d.patch1Start = getPatch1Start(d);
            d.patch2Start = getPatch2Start(d);
            d.patch1Stop = getPatch1Stop(d);
            d.patch2Stop = getPatch2Stop(d);
            d.patchSpeed = getPatchSpeed(d);
            d.tDur = getTDur(d);
            d.rewardVol = getRewardVol(d);
            d.conditionSwitch = getConditionSwitch(d);
            d.conditionSwitchList = getConditionSwitchList(d);
            d.SOA = getSOA(d);
            d.fixDuration = getFixDuration(d);

            if ~isempty(d.eye)
                d.eye_rm = rmBlinkSaccade(d);
                d.switchTime = getSwitchTime(d);
            end

            d.keyPressTime = getKeyPressTime(d);

        end

        function tDur = getTDur(d)
            tDur = d.meta.cic.tDur('time',Inf,'trial',1).data;
        end
        
        function fixDuration = getFixDuration(d)
            % fixDuration = d.meta.cic.fixDuration('time',Inf,'trial',1).data;
            fixDuration = d.meta.fixstim.fixDuration('time',Inf).data;
        end
        % function nRepPerCond = getNRepPerCond(d)
        %     nRepPerCond = d.meta.cic.nRepPerCond('time',Inf,'trial',1).data;
        % end
        
        function rewardVol = getRewardVol(d)
            rewardVol = d.meta.cic.rewardVol('time',Inf,'trial',1).data;
        end

        function conditionSwitchList = getConditionSwitchList(d)
            conditionSwitchList = d.meta.cic.conditionSwitch('time',Inf,'trial',1).data;
        end
        
        function conditionSwitch = getConditionSwitch(d)
            conditionSwitch = d.meta.patch1.conditionSwitch('time',Inf).data;
        end
        
        function SOA = getSOA(d)
            SOA = d.meta.cic.SOA('time',Inf,'trial',1).data;
        end

        function speed = getPatchSpeed(d)
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

        function patch1Dir = getPatch1Dir(d)
            patch1Dir = d.meta.patch1.direction('time',Inf).data;
        end

        function patch2Dir = getPatch2Dir(d)
            %patch2Dir = d.meta.patch2.direction('time',Inf).data; %NG in grating stim
            patch2Dir = d.patch1Dir+180*(1-d.congruent);
        end

        function patchDirList = getPatchDirList(d)
            patchDirList = unique(d.patch1Dir);
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

        function patch1Start = getPatch1Start(d)
            %Patch1Start = d.meta.patch1.on('time',Inf).data; %ms
            patch1Start = 1e3*(d.meta.patch1.startTime('time',Inf) - d.meta.cic.firstFrame('time',Inf))'; %ms
            patch1Start(patch1Start<0) = NaN;
            % .on is when you request neurostim to start a stimulus object
            % .startTime is when neurostim actually logs that the stimulus has started.
        end

        function patch2Start = getPatch2Start(d)
            %patch2Start = 1e3*(d.meta.patch2.startTime('time',Inf)- d.meta.cic.firstFrame('time',Inf))'; %ms
            patch2Start = 1e3*(d.meta.patch2.color('time',Inf).time - d.meta.cic.firstFrame('time',Inf))'; %ms
            patch2Start(patch2Start<0) = NaN;
        end

        function patch1Stop = getPatch1Stop(d)
            %patch1Stop = 1e3*(d.meta.patch1.stopTime('time',Inf)- d.meta.cic.firstFrame('time',Inf))'; %ms
            patch1Stop = 1e3*(d.meta.patch1.color('time',Inf).time - d.meta.cic.firstFrame('time',Inf))'; %ms
            patch1Stop(patch1Stop<0) = NaN; %ms
        end

        function patch2Stop = getPatch2Stop(d)
            patch2Stop = 1e3*(d.meta.patch2.stopTime('time',Inf)- d.meta.cic.firstFrame('time',Inf))'; %ms
            patch2Stop(patch2Stop<0) = NaN; %ms
        end

        function keyPressTime = getKeyPressTime(d)
            %get time of key press from the onset of each trial [ms]

            t0 =  d.meta.cic.firstFrame('time',Inf);

            [time,trial,frame,keyTmp] = d.meta.keypress.keyIx('time',Inf);
            key = cell2mat(keyTmp);
            ignoreTrial = isnan(key);
            keepInd = find(~ignoreTrial);
            time = time(~ignoreTrial);
            trial = trial(~ignoreTrial);
            frame = frame(~ignoreTrial);
            key = key(~ignoreTrial);

            keyPressTime = cell(d.numTrials, 1);
            for itr = 1:d.numTrials
                keyPressTime{itr} = 1e3*(time(trial == itr) - t0(itr));
                if isempty( keyPressTime{itr})
                     keyPressTime{itr} = NaN;
                end
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
            patch1On = d.meta.patch1.on('time',Inf).data; %intended time, not the actual time
            patch2On = d.meta.patch2.on('time',Inf).data; %intended time, not the actual time
            numAbort(1) = sum((1e3*trialEndTime <= patch1On) .* (d.complete == false)); %before 1st patch
            numAbort(2) = sum((1e3*trialEndTime > patch1On) .* (1e3*trialEndTime <= patch2On) .* (d.complete == false)); %between 1st and 2nd patches
            numAbort(3) = sum((1e3*trialEndTime > patch2On) .* (d.complete == false)); %after 2nd patch
            fprintf('# aborted trials before 1st patch: %d, between patches: %d, after 2nd patch: %d\n', numAbort(1), numAbort(2), numAbort(3));
        end

        function eye_rm = rmBlinkSaccade(d)
            for itr = 1:d.numTrials
                tmp = d.eye(itr).rmBlinks('dt',median(diff(d.eye(itr).t)), 'duration', 0.01, 'debug', false); %marmodata/+marmodata/@eye/rmBlinks.m
                eye_rm(itr,1) = tmp.rmSaccades('debug',false,'sargs',{'accthresh', d.accThresh}); %marmodata/+marmodata/@eye/rmSaccades.m
                %cf. fitKernel/selectSaccades
                %close all
            end
        end


        function [invalid, reason] = getInvalidBhvTrials(d)
            invalid = zeros(d.numTrials,3);
            reason = cell(d.numTrials,1);

            [~,~,vel_patchDir, avgGain_patchDir, avgAngdiff] = d.getSwitchTime;

            tmargin = d.fixDuration;

            for itr = 1:d.numTrials
                tidx = find((1e3*d.eye(itr).t > d.patch1Start(itr)-tmargin) & (1e3*d.eye(itr).t < d.patch2Stop(itr) + tmargin));%max(d.eye(itr).t)))
                if isempty(tidx); continue; end;

                t = 1e3*d.eye(itr).t(tidx); %ms
                t_r = t(1):1e3*1/d.fs_r:t(end);

                if sum(vel_patchDir{itr}(t_r < d.patch2Start(itr)) <0)>=1
                    invalid(itr,1) = 1;
                    reason{itr} = 'negative velocity during patch1';
                end

                if avgGain_patchDir(itr) < d.gainThresh
                    invalid(itr,2) = 1;
                    reason{itr} = [reason{itr} ', gain <' num2str(d.gainThresh)];
                end

                if mean(abs(avgAngdiff(itr,:))) > d.angdiffThresh
                    invalid(itr,3) = 1;
                    reason{itr} = [reason{itr} ', angluar diff > ' num2str(d.angdiffThresh)];
                end
            end

            disp([ num2str(sum(invalid(:,1))) ' negative velocity, '  num2str(sum(invalid(:,2))) ...
                ' too small gain, '   num2str(sum(invalid(:,3))) ' too large angular diff'])
        end


        function [switchTimes, switchTime1st, vel_patchDir, avgGain_patchDir, avgAngdiff] = getSwitchTime(d)
            switchTimes = cell(d.numTrials, 1); %ms %all detected movement switches during a trial
            switchTime1st = NaN(d.numTrials,1); %ms %earliest switch to 2nd patch
            vel_patchDir = cell(d.numTrials, 1); % deg/s
            avgGain_patchDir = nan(d.numTrials, 1);
            avgAngdiff = nan(d.numTrials, 2); %[deg]
            for itr = 1:d.numTrials

             tmargin = d.fixDuration;%ms
              tidx = find((1e3*d.eye(itr).t > d.patch1Start(itr) - tmargin) & (1e3*d.eye(itr).t < d.patch2Stop(itr) + tmargin));%max(d.eye(itr).t)));
            if isempty(tidx); continue; end;

                t = 1e3*d.eye(itr).t(tidx); %ms

                for ii = 1:2
                    switch ii
                        case 1
                            x = d.eye_rm(itr).x(tidx);
                            signal = x;
                        case 2
                            y = d.eye_rm(itr).y(tidx);
                            signal = y;
                    end

                    okIdx = find(~isnan(signal));
                    if isempty(okIdx)
                        continue;
                    end
                    signal = interp1(okIdx, signal(okIdx), 1:numel(signal), 'nearest','extrap')'; %extrapolate r to remove NaNs

                    signal_f = analysis.src.lowpassFilter(signal, 1e-3*t, d.cutoffFreq);

                    t_r = t(1):1e3*1/d.fs_r:t(end);
                    signal_rf = interp1(t, signal_f,t_r)'; %filtered and resampled
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
                vel_patchDir{itr} = sqrt(vx_rf.^2+vy_rf.^2) .* cos(theta - pi/180*d.patch1Dir(itr));
                patch1Tidx = find((t_r>d.patch1Start(itr)) & (t_r<d.patch2Start(itr)));
                patch2Tidx = find((t_r>d.patch2Start(itr)) & (t_r<d.patch2Stop(itr)));
                avgAngdiff(itr, 1) = median(180/pi*analysis.src.angdiff(theta(patch1Tidx), repmat(pi/180*d.patch1Dir(itr), [numel(patch1Tidx),1])));
                avgAngdiff(itr, 2) = median(180/pi*analysis.src.angdiff(theta(patch2Tidx), repmat(pi/180*d.patch2Dir(itr), [numel(patch2Tidx),1])));
                avgGain_patchDir(itr) = median(abs(vel_patchDir{itr}(t_r>d.patch1Start(itr))))/d.patchSpeed;

                %% detection of switch in eye movement direction
                [pks, peakTidx_c] = findpeaks(-abs(vel_patchDir{itr}), 'MinPeakProminence',median(abs(vel_patchDir{itr}(t_r>d.patch1Start(itr)))));%d.prominenceThresh);
                peakTidx = [];
                for jj = 1:numel(pks)
                    if pks(jj) < 0.1 && peakTidx_c(jj) > d.dt_ba*d.fs_r+1 && peakTidx_c(jj) < numel(t_r)-d.dt_ba*d.fs_r && ...
                            vel_patchDir{itr}(peakTidx_c(jj)-d.dt_ba*d.fs_r)*vel_patchDir{itr}(peakTidx_c(jj)+d.dt_ba*d.fs_r)<0
                        peakTidx = [peakTidx peakTidx_c(jj)];
                    end
                end
                if ~isempty(peakTidx)
                    switchTimes{itr} = t_r(peakTidx); %all switches
                else
                    switchTimes{itr}= NaN;
                end

                %1st switch to 2nd patch from positive to negative velocity
                if ~isempty(peakTidx)
                    thisIdx = find( (vel_patchDir{itr}(peakTidx-d.dt_ba*d.fs_r)'>0) .* (vel_patchDir{itr}(peakTidx+d.dt_ba*d.fs_r)'<0) ...
                        .* (t_r(peakTidx)>d.patch2Start(itr)) .* (t_r(peakTidx)<d.patch2Stop(itr)));
                    if ~isempty(thisIdx)
                        switchTime1st(itr) = t_r(min(peakTidx(thisIdx))); 
                    end
                end
            end
        end

        function eyeKeyConsistency(d)
            keySwitchTime = d.keyPressTime;
            [~, switchTime1st] = d.getSwitchTime;

            keySwitched =  cellfun(@(x)~isnan(x), keySwitchTime);
            eyeSwitched = ~isnan(switchTime1st);

            theseTrials = find(d.complete.*~d.congruent);
            both = keySwitched(theseTrials) .* eyeSwitched(theseTrials);
            eyeOnly = (keySwitched(theseTrials)==0) .* (eyeSwitched(theseTrials) == 1);
            keyOnly = (keySwitched(theseTrials)==1) .* (eyeSwitched(theseTrials) == 0);
            consistent =  both + (keySwitched(theseTrials)==0) .* (eyeSwitched(theseTrials) == 0); %Hasse and Tsao 2020 eLife

            fprintf('eye and key: %d/%d trials (%2.2f%%) \nkey only: %d/%d trials (%2.2f%%) \neye only: %d/%d trials (%2.2f%%) \nconsistent: %d/%d trials (%2.2f%%) \n', ...
                sum(both), numel(theseTrials), 100*sum(both)/numel(theseTrials), ...
                sum(keyOnly), numel(theseTrials), 100*sum(keyOnly)/numel(theseTrials), ...
                sum(eyeOnly), numel(theseTrials), 100*sum(eyeOnly)/numel(theseTrials), ...
                sum(consistent), numel(theseTrials), 100*sum(consistent)/numel(theseTrials));
        end

        %% behavioural stats
        function f = statEyeSwitch(d)
            [~, switchTime1st] = d.getSwitchTime;
           
            %switch rate
            disp('#trials with eye switch after 2nd patch:');
           for icond = d.conditionSwitchList
               switch icond
                   case 0
                       condName = 'flash suppression';
                   case 1
                       condName = 'physical alteration';
                   case 2
                       condName = 'congruent';
               end
               trials = find(d.conditionSwitch == icond & d.complete);
               numAllTrials = numel(trials);
               numEyeSwitchTrials = sum(~isnan(switchTime1st(trials)));
               fprintf('%s: %d/%d trials (%.2f%%)\n', condName, numEyeSwitchTrials, numAllTrials, 100*numEyeSwitchTrials/numAllTrials);
           end


           %latency
           latencyFS = switchTime1st(d.conditionSwitch==0.*d.complete) - d.patch2Start(d.conditionSwitch==0.*d.complete);
           latencyPA = switchTime1st(d.conditionSwitch==1.*d.complete) - d.patch2Start(d.conditionSwitch==1.*d.complete);
           fprintf('mean eye switch latency across trials in FS: %.1f, PA: %.1f [ms]\n', nanmean(latencyFS), nanmean(latencyPA));

           f=figure;
           histogram(latencyFS, 20);
           hold on
           histogram(latencyPA, 20);
           grid on
           vline(nanmean(latencyFS),gca,[],'b')
           vline(nanmean(latencyPA),gca,[],'r')
           legend('FS','PA')
           xlabel('eye switch latency [ms]')
           ylabel('#trials')

        end

        %% functions for visualization
        function fig = plotSwitchByPatchDir(d)
            [~, switchTime1st] = d.getSwitchTime;
            % switched = cellfun(@(x)~isnan(x), switchTimes);
            switched = ~isnan(switchTime1st);

            nCompleteTrials = zeros(numel(d.patchDirList), 2, 2);
            nSwitchedTrials = zeros(numel(d.patchDirList), 2, 2);
            for idir = 1:numel(d.patchDirList)
                for ifirstRed = 0:1
                    for icongruent = 0:1
                        completeTrials = (d.patch1Dir == d.patchDirList(idir)) .* ...
                            (d.redFirst==ifirstRed) .* (d.complete==true) .*  (d.congruent == icongruent);
                        nCompleteTrials(idir, ifirstRed+1, icongruent+1) = sum(completeTrials); %supposed to be all d.numTrials/d.numConds

                        nSwitchedTrials(idir, ifirstRed+1, icongruent+1) = sum((switched==1) .* completeTrials);
                    end
                end
            end
            fig = figure;
            ax(1)=subplot(211);
            plot(d.patchDirList, nSwitchedTrials(:,1,1)./nCompleteTrials(:,1,1), '-*','Color','b'); hold on
            plot(d.patchDirList, nSwitchedTrials(:,2,1)./nCompleteTrials(:,2,1), '-o','Color','r');
            title('incongruent');
            axis padded;

            ax(2)=subplot(212);
            plot(d.patchDirList, nSwitchedTrials(:,1,2)./nCompleteTrials(:,1,2), '-*','Color','b'); hold on
            plot(d.patchDirList, nSwitchedTrials(:,2,2)./nCompleteTrials(:,2,2), '-o','Color','r');
            title('congruent'); 
          
            axis padded;
            linkaxes(ax(:));
            legend('blueFirst','redFirst');
            xlabel('stimulus direction [deg]'); ylabel('eye switch rate');
        end

        function fig = plotSingleTrialEye(d,itr)
            %summary plot for eye movement per trial, showing
            %raw eye position
            %eye position after removal of blinks and saccades
            % eye velocity along the 1st patch direction

            winSize = 5;
            
            [switchTimes, switchTime1st, vel_patchDir] = d.getSwitchTime;
            
           tmargin = d.fixDuration;%ms
            tidx = find((1e3*d.eye(itr).t > d.patch1Start(itr)-tmargin) & (1e3*d.eye(itr).t < d.patch2Stop(itr) + tmargin));%max(d.eye(itr).t)));
            if isempty(tidx); return; end

            t = 1e3*d.eye(itr).t(tidx); %s
            t_r = t(1):1e3/d.fs_r:t(end);

            [~, Patch1Starttidx] = min(abs(t - d.patch1Start(itr)));
            [~, Patch2Starttidx] = min(abs(t - d.patch2Start(itr)));

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
            vline(d.patch2Start(itr),gca,'-',thisColor);
            vbox(d.patch1Start(itr),d.patch1Stop(itr),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.patch2Start(itr),d.patch2Stop(itr),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            vline(switchTimes{itr}, gca, '--','g'); vline(switchTime1st(itr), gca, '-','g');
            vline(d.keyPressTime{itr}, gca, '-','y');
            title(sprintf('tr:%d', itr))
            xlabel('time [ms]'); ylabel('raw eye position [deg]')
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(222);
            plot(t, x_rm, t, y_rm);
            vline(d.patch2Start(itr),gca,'-',thisColor);
            vbox(d.patch1Start(itr),d.patch1Stop(itr),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.patch2Start(itr),d.patch2Stop(itr),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            vline(switchTimes{itr}, gca, '--','g'); vline(switchTime1st(itr), gca, '-','g');
            vline(d.keyPressTime{itr}, gca, '-','y');
            xlabel('time [ms]'); ylabel('after removal of saccade [deg]');
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(223);
            %plot(1e3*t, r, 'k',1e3*t, r_f, 'g');
            plot(t_r, vel_patchDir{itr});
            vline(d.patch2Start(itr),gca,'-',thisColor);
            vbox(d.patch1Start(itr),d.patch1Stop(itr),gca,[d.redFirst(itr) 0.0 1-d.redFirst(itr) .1])
            vbox(d.patch2Start(itr),d.patch2Stop(itr),gca,[1-d.redFirst(itr) 0.0 d.redFirst(itr) .1])
            hline(0);
            vline(switchTimes{itr}, gca, '--','g'); vline(switchTime1st(itr), gca, '-','g');
            vline(d.keyPressTime{itr}, gca, '-','y');
            xlabel('time [ms]');
            ylabel('velocity along 1st patch');
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end

            subplot(224)
            quiver(0,0, winSize*cos(d.patch1Dir(itr)*pi/180), winSize*sin(d.patch1Dir(itr)*pi/180),'linewidth',1); hold on;
            try
                x0 = x_rm(Patch1Starttidx);
                y0 = y_rm(Patch1Starttidx);
                [~,switchTime1stIdx] = min(abs(t - switchTime1st(itr)));
                scatter(x_rm(Patch1Starttidx:end) - x0, y_rm(Patch1Starttidx:end) - y0, 1, t(Patch1Starttidx:end));
                plot(x_rm(Patch2Starttidx)- x0, y_rm(Patch2Starttidx)- y0, 'o', 'Color',thisColor);
                if ~isnan(switchTime1st(itr))
                    plot(x_rm(switchTime1stIdx)- x0, y_rm(switchTime1stIdx)- y0, 'go');
                end
            end
            colormap("cool");
            xlim([-winSize winSize]);  ylim([-winSize winSize]);
            vline(0,gca,'--',[.5 .5 .5 .5]); hline(0,gca,'--',[.5 .5 .5 .5]);
            xlabel('x [deg]'); ylabel('y [deg]');
            axis square;
            if d.complete(itr)==0; set(gca,'xcolor','r','ycolor','r');end
        end

        function fd = checkDroppedFrames(d)
            %% under construction

            fd.trial=d.meta.cic.frameDrop.trial';
            fd.time=d.meta.cic.frameDrop.time';
            startTime = d.meta.cic.firstFrame('time',inf).time;

            dt = [];
            tr = unique(fd.trial);
            figure('Name','Dropped frames')
            subplot(3,1,1:2)
            for a = tr'
                x = fd.time(fd.trial==a) - startTime(a); % times for this trial
                dt = [dt; diff(x)];
                plot(x,a*ones(size(x)),'ko')
                hold on
            end
            %xlim([0 max(fd.trialTime)])
            ylabel('Trial')
            xlabel('trialTime (ms)')

            subplot(3,1,3)
            histogram(dt,0:10:2000)
            xlabel('delta (ms)')
            ylabel('count')

            %% stats on FS and congruent trials
            theseTrials = find(d.complete.*(d.conditionSwitch~=1));
            patch2Dur = d.patch2Stop(theseTrials) - d.patch2Start(theseTrials);
            disp(['median duration of patch2 in FS/congruent: ' num2str(median(patch2Dur))]);
            disp(['median dropped frames in FS/congruent: ' num2str(sum(ismember(fd.trial, theseTrials))/numel(theseTrials))]);

        end

    end
end