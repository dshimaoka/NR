
deg = 0:15:180;

%% test1: change directions of both eyes
% nystagmus_directions_twoEyes('twoEyes','phaseSpeed',10, 'ori1List',deg, ...
%     'nRepPerCond',1, 'sigma',2, 'frequency',.5,'tDur',1e4);

%% test2: change difference between eyes (left eye red) is fixed to 0)
% nystagmus_directions_oneEye('oneEye','phaseSpeed',10, 'ori1List',deg, ...
%     'nRepPerCond',1, 'sigma',2, 'frequency',.5,'tDur',1e4);



%% analysis
params.paths = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\recording\nystagmus_20220820';%'/home/marmolab/data/2022/08/20/';
params.files = 'twoEyes.pursuit2D.113116.mat';
%params.files = 'oneEye.pursuit2D.113617.mat';

d = marmodata.mdbase('path',params.paths,'file',params.files,'loadArgs',{'loadEye',true});

%% show raw traces in x & y of every conditions
load(fullfile(params.paths, params.files), 'c');
stimDir = get(c.gabor1.prms.orientation,'atTrialTime',inf);
successTr = ~isnan(get(c.choice.prms.keyIx,'atTrialTime',inf));

figure('position',[0 0 1920 1080]);
for itr = 1:d.numTrials
    if ~successTr(itr)
        continue;
    end
    ax(itr)=subplot(3,5,find(deg == stimDir(itr)));
    tidx = find((d.eye(itr).t > 0) & (d.eye(itr).t < get(c.gabor1.prms.duration,'atTrialTime',inf,'trial',itr)*1e-3)); 
    plot(d.eye(itr).t(tidx), d.eye(itr).x(tidx), d.eye(itr).t(tidx), d.eye(itr).y(tidx));
    title([num2str(stimDir(itr)) ' deg']);
end
linkaxes(ax(:));
xlabel('time [s]')
ylabel('eye position [deg]')
legend('x','y');
saveas(gcf,fullfile(params.paths,['1d_' params.files(1:end-4) '.png']));

%% remove saccades then show traces in x-y 2D space
ey_rmSaccades = [];
for itr = 1:d.numTrials
    eye_rmSaccades(itr) = d.eye(itr).rmSaccades;
    close
end

figure('position',[0 0 1920 1080]);
for itr = 1:d.numTrials
    if ~successTr(itr)
        continue;
    end
    ax(itr)=subplot(3,5,find(deg == stimDir(itr)));
    tidx = find((eye_rmSaccades(itr).t > 0) & (eye_rmSaccades(itr).t < get(c.gabor1.prms.duration,'atTrialTime',inf,'trial',itr)*1e-3)); 
    plot(eye_rmSaccades(itr).t(tidx), eye_rmSaccades(itr).x(tidx), eye_rmSaccades(itr).t(tidx), eye_rmSaccades(itr).y(tidx));
    title([num2str(stimDir(itr)) ' deg']);
end
linkaxes(ax(:));
xlabel('time [s]')
ylabel('eye position [deg]')
legend('x','y');
saveas(gcf,fullfile(params.paths,['1d_' params.files(1:end-4) '_rmSaccades.png']));


%% figure2: x-y 2D space
figure('position',[0 0 1920 1080]);
for itr = 1:d.numTrials
    if ~successTr(itr)
        continue;
    end
    ax(itr)=subplot(3,5,find(deg == stimDir(itr)));
    tidx = find((eye_rmSaccades(itr).t > 0) & (eye_rmSaccades(itr).t < get(c.gabor1.prms.duration,'atTrialTime',inf,'trial',itr)*1e-3)); 
    plot(eye_rmSaccades(itr).x(tidx) - eye_rmSaccades(itr).x(tidx(1)), ...
        eye_rmSaccades(itr).y(tidx) - eye_rmSaccades(itr).y(tidx(1)));
    hold on
    quiver(0,0,-10*sin(pi/180*stimDir(itr)), 10*cos(pi/180*stimDir(itr)));
    title([num2str(stimDir(itr)) ' deg']);
end
xlabel('x after removal of saccade [deg]')
ylabel('x after removal of saccade [deg]')
linkaxes(ax(:));
axis equal
legend('eye','stimulus')
saveas(gcf,fullfile(params.paths,['1d_' params.files(1:end-4) '_rmSaccades_xy.png']));
