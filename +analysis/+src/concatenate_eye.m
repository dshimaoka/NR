function [eyeData_cat, meta_cat] = ...
    concatenate_eye(eyeData, dd)
%[eyeData_cat] = concatenate_eye(dd)
% returns mamordata.eye concatenated across trials

%[eyeData_cat, meta_cat] = concatenate_eye(dd)
% additionally reutns saccade and blink times computed in eyelink

nTrials = length(eyeData);


edfnames = {'STARTSACC','ENDSACC','STARTBLINK','ENDBLINK'};
c = cell(length(edfnames),1);
meta_cat = cell2struct(c,edfnames);

%% concatenate across trials as in my Cell Rep paper 2018
%does parea=0 mean blink?
%is there temporal gap between trials?
parea_cat = [];
pwdth_cat = [];
phght_cat = [];
x_cat = [];
y_cat = [];
t_cat = [];
%dd_cat = marmodata.cuesaccade;


for itr = 1:nTrials
    x_cat = cat(1, x_cat, eyeData(itr).x);
    y_cat = cat(1, y_cat, eyeData(itr).y);
    parea_cat = cat(1, parea_cat, eyeData(itr).parea);
    phght_cat = cat(1, phght_cat, eyeData(itr).phght);
    pwdth_cat = cat(1, pwdth_cat, eyeData(itr).pwdth);

    if isempty(t_cat)
        t0 = eyeData(itr).t(1);
    elseif ~isempty(eyeData(itr).t)
        t0 =  max(t_cat)-eyeData(itr).t(1)+eyeData(itr).dt;
    end
    if ~isempty(eyeData(itr).t)
        t_cat = cat(1, t_cat, eyeData(itr).t+t0);
    end

    for iedf = 1:length(edfnames)
        meta_cat.(edfnames{iedf}) = cat(1, meta_cat.(edfnames{iedf}), ...
            (dd.meta.edf.(edfnames{iedf})('trial',itr).time)' ...
            - dd.meta.cic.firstFrame('trial',itr).time+t0);
    end
end
t_cat = t_cat(~isnan(t_cat));
[t_cat, ix] = unique(t_cat);
eyeData_cat = marmodata.eye(t_cat, x_cat(ix), y_cat(ix), pwdth_cat(ix), phght_cat(ix));


%% check STARTBLINK/ENDBLINK
if length(meta_cat.STARTBLINK) >  length(meta_cat.ENDBLINK)
    meta_cat.ENDBLINK = [meta_cat.ENDBLINK; eyeData_cat.t(end)];
elseif length(meta_cat.STARTBLINK) <  length(meta_cat.ENDBLINK)
    meta_cat.STARTBLINK = [meta_cat.STARTBLINK; eyeData_cat.t(1)];
end


if ~isempty(find(meta_cat.ENDBLINK - meta_cat.STARTBLINK < 0))
    error('ENDBLINK comes earlier than STARTBLINK');
end

%% check STARTSACC/ENDSACC 28/1/22
if length(meta_cat.STARTSACC) ~= length(meta_cat.ENDSACC)%41
    nSaccs = min(length(meta_cat.STARTSACC), length(meta_cat.ENDSACC));
    ngIdx=find(meta_cat.ENDSACC(1:nSaccs)-meta_cat.STARTSACC(1:nSaccs)<0);

    if isempty(ngIdx)
        meta_cat.STARTSACC = meta_cat.STARTSACC(1:nSaccs);
        meta_cat.ENDSACC = meta_cat.ENDSACC(1:nSaccs);
        %else
        % FILL ME?
    end
end
idx = find(meta_cat.STARTBLINK-meta_cat.ENDBLINK~=0);
meta_cat.STARTBLINK = meta_cat.STARTBLINK(idx);
meta_cat.ENDBLINK = meta_cat.ENDBLINK(idx);

assert(isempty(find(meta_cat.ENDSACC-meta_cat.STARTSACC<0)));
okSacc = find(meta_cat.ENDSACC-meta_cat.STARTSACC>0);
meta_cat.STARTSACC = meta_cat.STARTSACC(okSacc);
meta_cat.ENDSACC = meta_cat.ENDSACC(okSacc);
