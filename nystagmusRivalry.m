function nystagmusRivalry(subject,varargin)
% created from pursuit2D.m  https://github.com/nPrice-lab/pursuit
%
% mouse as a eye position
% mouse position is registered and fixed by a left click
% it stays at the same position until further clicking
%
%
%
%
% After paradigm starts:
% - Press Enter to see your eye - adjust focus of camera
% - Press 'a' for automatic luminance adjustment, then `c` to calibrate.
% - Press 'space' when you first fixate and 'Enter' when
%   you've finished looking at the 9 spots.
% - Press Esc (possibly twice) when you're ready to start protocol.
% - Press 'a' or 'z' to start a new trial
%   *********** Press "Esc" twice to exit a running experiment ************
%
%
% Trial structure:
% - acquire and hold fixation for a fixed (300ms) duration
% - presentation of a first patch, followed by a second patch after SOA
% - maitaining fixation within radius gets rewarded
%

%% below from pursuit.m
% Implementation
% - initial static fixation spot is one "stimulus"
% - moving spot is drawn by a second "stimulus", with position controlled by
% a third stimulus function - targXY.m - whose only role is to update the
% XY position.
% - movement starts at pre-defined time after fixation is acquired
%
%
% Also
% - restart Matlab before you collect data in a new session.
% - close all potentially demanding programs like Firefox, Slack, GitKraken etc.
%       File manager and terminals are fine.
% - turn off EyeLink monitor or push it to edge of table (after you know it's working)
% - close door to hallway
% - have lights in anteroom on, in testing room off
%
% Auto backup recent data:
% - Ctrl-Alt-T (pop a terminal)
% - syncit (mount syncitium) It may ask "password for nic", which is the computer password - (Hint: monkey)
%   it will also ask "Enter syncitium username", which is your standard Monash username and then password.
% - rsyncData (push data)
%
% Check frame drops:
% fName = '/home/marmolab/data/2022/11/29/tst.pursuit2D.220253.mat';
% load(fName)
% fd=get(c.prms.frameDrop,'struct',true);
%
%
%

%% PARAMETER DEFINITIONS

if ~exist('subject','var')
    error('No subject name provided. Type ''help facecal'' for usage information.');
end

validateattributes(subject,{'char'},{'nonempty'},'','subject',1);

% parse arguments...
p = inputParser();
p.KeepUnmatched = true;
p.addRequired('subject',@(x) validateattributes(x,{'char'},{'nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('tDur',4000,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % trial duration from onset of first patch (ms)
p.addParameter('nRepPerCond',10,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % number of repeats of each condition
p.addParameter('rewardVol',0.035,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); % adopted from OcuFol
p.addParameter('conditionSwitch', [0 1 2], @(x) validateattributes(x,{'numeric'},{'vector','nonempty'}));
%conditionSwitch = 0: binocular flash suppression
%conditionSwitch = 1: physical alteration
%conditionSwitch = 2: congruent direction between two patches

%patch stimuli
p.addParameter('patchType','rdp',@(x) validateattributes(x,{'char'},{'nonempty'})); %rdp or grating
p.addParameter('dirList_first',[0], @(x) validateattributes(x,{'numeric'},{'vector','nonempty'})); %direction(s) of the first patch [deg] 0: left to right, 90: bottom to top
p.addParameter('speed',12, @(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); %[(visual angle in deg)/s]
p.addParameter('radius',5, @(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); %aperture size [deg]
p.addParameter('SOARange', [1000 1001], @(x) validateattributes(x,{'numeric'},{'vector','nonempty'})); %stimulus onset after the end of fixation

p.parse(subject,varargin{:});
args = p.Results;

%% fixed parameters
fixDuration = 300; % [ms] minimum duration of fixation to initiate patch stimuli
fixationDeadline = 5000; %[ms] maximum time to initiate a trial 
iti = 1000; %[ms] inter trial interval

%luminance correction
% redLuminance = 171/255; %Fraser ... Miller 2023
redLuminance = 0.33; %DS office

%RDP
dotSize = 1;%5; %dot size [pix]
nrDots = 200; %number of dots

%grating
frequency = 0.5; %spatial frequency in cycles per visual angle in degree (not pixel) %Kapoor 2022

import neurostim.*
commandwindow;

% total trial number
% args.nRepPerCond * numel(args.dirList_first) * 2 * 2 % direction x (congruent / incongruent) * (red/blue)


%% ========= Specify rig configuration  =========

%Create a Command and Intelligence Centre object (the central controller for everything). Here a cic is returned with some default settings for this computer, if it is recognized.
c = marmolab.rigcfg('debug',args.debug, p.Unmatched); % set to false to save githash at start of each experiment!
%c = myRig;
c.paradigm = 'nystagmusRivalry';
c.addProperty('fixDuration', fixDuration);
c.addProperty('jitteredSOA',[]);
c.jitteredSOA = plugins.jitter(c,{args.SOARange(1), args.SOARange(2)}); 
c.addProperty('tDur',args.tDur);
c.screen.color.background = 0.5+[0 0 0];
c.addProperty('redLuminance', redLuminance);

if ~args.debug % log git hash
    hash = marmolab.getGitHash(fileparts(mfilename('fullpath')));
    c.githash('NR.git') = hash;
end

%Make sure there is an eye tracker (or at least a virtual one)
if isempty(c.pluginsByClass('eyetracker'))
    e = neurostim.plugins.eyetracker(c);      %Eye tracker plugin not yet added, so use the virtual one. Mouse is used to control gaze position (click)
    e.useMouse = true;
end


%% ============== Add stimuli ==================

%% Fixation dot
f = stimuli.fixation(c,'fixstim');    % Add a fixation stimulus object (named "fix") to the cic. It is born with default values for all parameters.
f.shape = 'CIRC';               %The seemingly local variable "f" is actually a handle to the stimulus in CIC, so can alter the internal stimulus by modifying "f".
f.size = 0.25; % units?
f.color = [1 1 1];
f.on=0;                         % What time should the stimulus come on? (all times are in ms)
f.duration = '@fixbhv.startTime.fixating+cic.fixDuration'; % Show spot briefly after fixation acquired
f.X = 0;
f.Y = 0;

%% RDP
s = RandStream('mt19937ar');

nrConds = 2;
fm = cell(nrConds,1);
for ii = 1:nrConds
    
    reset(s, 1);%args.rngSeed);

    stimName = ['patch' num2str(ii)];
    %patch1: presented 1st, patch2: presented 2nd after SOA
    
    if strcmp(args.patchType, 'rdp')
        fm{ii} = neurostim.stimuli.rdp(c,stimName);
        %rdp specific parameters
        fm{ii}.maxRadius =  args.radius;%maximum radius of aperture (px)
        fm{ii}.speed = args.speed; %dot speed (deg
        fm{ii}.size = dotSize; %dot size [px]
        fm{ii}.type = 0; %square dot
        fm{ii}.nrDots = nrDots;
        fm{ii}.coherence = 1; %dot coherence [0-1]
        fm{ii}.motionMode = 1; %linear
        fm{ii}.lifetime = 50;%lifetime of dots (in frames)
        fm{ii}.dwellTime = 1;
        fm{ii}.coordSystem = 0; %polar coordinates
        fm{ii}.noiseMode = 0; %proportion
        fm{ii}.noiseDist = 1; %uniform
        fm{ii}.direction = 0;
    elseif strcmp(args.patchType,'grating')
        fm{ii} = neurostim.stimuli.gabor(c, stimName);
        fm{ii}.width = 2*max(args.radius);
        fm{ii}.height = fm{ii}.width;
        fm{ii}.sigma = args.radius;
        fm{ii}.mask = 'CIRCLE';
        fm{ii}.frequency = frequency;
        fm{ii}.contrast = 1;
        fm{ii}.flickerMode = 'sinecontrast';%'none'; %none makes the phase difference between patches more apparent
        fm{ii}.flickerFrequency = 0;
        fm{ii}.orientation = 0;
        fm{ii}.addProperty('direction',0);
        fm{ii}.addProperty('directionPolarity',0);
        fm{ii}.addProperty('speed',args.speed);  
        %TODO: align temporal phase between patches?
    end
    
    %common parameters across stim
    fm{ii}.X = 0;
    fm{ii}.Y = 0;
    fm{ii}.addProperty('frameRate', c.screen.frameRate);

end
fm{1}.addProperty('conditionSwitch', 1);

fm{1}.addProperty('redFirst',0);
%fm{1}.redFirst = plugins.jitter(c,{0, 1},'distribution','1ofN'); %NG always return 1
fm{1}.color = '@[cic.redLuminance*patch1.redFirst 0.0 1-patch1.redFirst 1]'; 
fm{2}.color = '@[1-patch1.redFirst 0.0 patch1.redFirst 1]'; 
fm{1}.on = '@fixstim.off'; %first stimulus
fm{2}.on = '@patch1.on + cic.jitteredSOA'; %2nd stimulus
fm{1}.duration = '@cic.tDur  - patch2.physicalAlteration * (cic.tDur - cic.jitteredSOA)';  %if physicalAlteration=1, terminate after jitteredSOA
fm{2}.duration = '@cic.tDur - cic.jitteredSOA';
fm{2}.addProperty('congruent', '@fix(patch1.conditionSwitch/2)'); %whether the second patch moves the same direction with the 1st patch
fm{2}.addProperty('physicalAlteration','@rem(patch1.conditionSwitch, 2)')

if strcmp(args.patchType, 'rdp')
    fm{2}.direction = '@patch1.direction+180*(1-patch2.congruent)';
elseif strcmp(args.patchType,'grating')
    fm{1}.orientation = '@mod(patch1.direction, 180) - 90';
    fm{2}.orientation = '@patch1.orientation';
    fm{1}.directionPolarity = '@-2*fix(patch1.direction/180) + 1';
    fm{2}.directionPolarity = '@(2*patch2.congruent-1) * patch1.directionPolarity';
    fm{1}.phaseSpeed = '@360*patch1.directionPolarity * patch1.speed * patch1.frequency /patch1.frameRate'; %[deg/frame]
    fm{2}.phaseSpeed = '@360*patch2.directionPolarity * patch2.speed * patch2.frequency /patch2.frameRate'; %[deg/frame] 
    %fm{2}.phase = '@patch1.phase + patch1.phaseSpeed*(patch1.frameRate+10*(1-patch2.physicalAlteration-patch2.congruent))*patch1.duration/1000 + 270*patch2.congruent'; %[deg] %works 1&2 not 0
    %fm{2}.phase = '@patch1.phase + patch1.phaseSpeed*patch1.frameRate*patch1.duration/1000 + 270'; %works in condSwitch=0(&2) not 1
    %fm{2}.phase = '@patch1.spatialPhase'; %NG0,1,2
    fm{1}.phase = '@mod(-patch1.phaseSpeed*patch1.frameRate*patch1.duration/1000 - 270 - 90*patch2.physicalAlteration, 360)';%works in condSwitch=0(&2) not 1
    fm{2}.phase = 0;  

    %% annulus mask ... under development
    % mm = neurostim.stimuli.convPoly(c, 'maskGrating');
    % mm.radius = args.radius;
    % mm.nSides = 20;
    % mm.filled = true;
    % %mm.linewidth = 5;
    % mm.color = [0 0 0 0];

    % mmo = neurostim.stimuli.convPoly(c, 'maskGrating_outer');
    % mmo.radius = args.radius+5;
    % mmo.nSides = 20;
    % mmo.filled = true;
    % %mm.linewidth = 5;
    % mmo.color = [1 1 1 0];

    % mmo = neurostim.stimuli.noiseradialgrid(c, 'maskGrating_outer');
    % mmo.nRadii = 2;
    % mmo.nWedges = 2;
    % mmo.parms = [.99 1];

end

%% ========== Add required behaviours =========
%Subject's 2AFC response to control inter-trial interval ... not necessary?
k = behaviors.keyResponse(c,'keypress');
k.from = '@patch2.off'; % end of patch
k.maximumRT= Inf;                   %Allow inf time for a response
k.keys = {'a'};%,'z'};
k.required = false; %   setting false means that even if this behavior is not successful (i.e. the wrong answer is given), the trial will not be repeated.

%Maintain gaze on the fixation point until the trial end
g = behaviors.fixate(c,'fixbhv');
g.from = fixationDeadline; % If fixation has not started at this time, move to the next trial
g.to = '@patch2.off'; %'@traj.off'; % stop tracking when trajectory ends
g.X = '@patch1.X'; %'@traj.X';
g.Y = '@patch1.Y'; %'@traj.Y';
g.tolerance = args.radius; % (deg) allowed eye position error - should be aiming to get this as small as possible
g.required = true; % This is a required behavior. Any trial in which fixation is not maintained throughout will be retried. (See myDesign.retry below)
g.failEndsTrial = true;
g.successEndsTrial = true; %cf. false in OcuFol


%% Turn off logging
stopLog(c.fixstim.prms.X);
stopLog(c.fixstim.prms.Y);
stopLog(c.patch1.prms.X);
stopLog(c.patch1.prms.Y);
stopLog(c.patch2.prms.X);
stopLog(c.patch2.prms.Y);
stopLog(c.fixbhv.prms.X);
stopLog(c.fixbhv.prms.Y);


%% ========== Specify feedback/rewards =========
% % Play a correct/incorrect sound for the 2AFC task
% plugins.sound(c);           %Use the sound plugin
% 
% % Add correct/incorrect feedback
% s= plugins.soundFeedback(c,'soundFeedback');
% s.add('waveform','correct.wav','when','afterTrial','criterion','@choice.correct');
% s.add('waveform','incorrect.wav','when','afterTrial','criterion','@ ~choice.correct');

if true && ~isempty(c.pluginsByClass('newera'))
  % add liquid reward... newera syringe pump  
  c.newera.add('volume',args.rewardVol,'when','AFTERTRIAL','criterion','@fixbhv.isSuccess');
end


%% Experimental design
c.trialDuration = Inf; %'@choice.stopTime';       %End the trial as soon as the 2AFC response is made.
% c.trialDuration = '@choice.stopTime + faces.duration'; %cuesaccde
% c.trialDuration = '@tarBr.startTime.fixating + tarBr.ps + 300'; %OcuFol
c.iti = iti;

%  Specify experimental conditions
% For threshold estimation, we'd just vary speed
myDesign=design('myFac');                      %Type "help neurostim/design" for more options.
facOutList = {'direction'}; % frequency = spatial frequency
facInList = {'dirList_first'};
for a = 1:length(facInList)
    myDesign.(sprintf('fac%d',a)).patch1.(facOutList{a}) = args.(facInList{a});
end

myDesign.fac2.patch1.redFirst = 1;%[0 1]; %whether to start with red or blue
myDesign.fac3.patch1.conditionSwitch = args.conditionSwitch;

myDesign.retry = 'RANDOM'; %'IMMEDIATE' or 'IGNORE';
myDesign.maxRetry = 4;%10;  % Each condition will be retried up to this many times.

a=1;
myBlk{a} = block('myBlock',myDesign);
myBlk{a}.nrRepeats = args.nRepPerCond; %params.nRepPerCond; %nRepeatsPerBlock;

c.eye.doTrackerSetupEachBlock = true; %KY disabled

%% Run the experiment.
% c.eye.clbMatrix = marmolab.loadCal(args.subject); %KY
% c.setPluginOrder('mov','blank','fix','tar','fWindow','sWindow'); %KY

c.subject = args.subject; %params.subj; %'NP';
c.run(myBlk{1}); %cf. KY c.run(myBlk,'nrRepeats',500);


% Quick check of framedrops
%
% load last file
% fd=get(c.prms.frameDrop,'struct',true);
% -OR-
% d = marmodata.mdbase('file',f,'loadArgs',{'loadEye',true});
% for a = 1:max(d.meta.cic.trial.data)
%   frDrop{a} = d.meta.cic.frameDrop('trial',a).data;
% end
% nDrop = cellfun('length',frDrop)
%
%% Check frame drops:
% fName = '/home/marmolab/data/2022/12/13/easyD.test.161617.mat'; %'/home/marmolab/data/2022/12/13/tst.pursuit2D.161358.mat'; %'/home/marmolab/data/2022/12/13/tst.pursuit2D.152612.mat' %'/home/marmolab/data/2022/11/29/tst.pursuit2D.220253.mat';
% load(fName)
fd=get(c.prms.frameDrop,'struct',true);

dt = [];
tr = unique(fd.trial);
figure('Name','Dropped frames')
subplot(3,1,1:2)
for a = tr'
    x = fd.trialTime(fd.trial==a); % times for this trial
    dt = [dt; diff(x)];
    plot(x,a*ones(size(x)),'ko')
    hold on
end
xlim([0 max(fd.trialTime)])
ylabel('Trial')
xlabel('trialTime (ms)')

subplot(3,1,3)
histogram(dt,0:10:2000)
xlabel('delta (ms)')
ylabel('count')