function signal_f = lowpassFilter(signal, t, cutoffFreq,order)
if nargin < 4
    order = 3;
end
if nargin < 3
    cutoffFreq = 1;%hz
end
fs = 1/median(diff(t));
Wn = cutoffFreq/(fs/2);
[b,a]=butter(order, Wn, 'low');

signal_c = filtfilt(b,a,cat(1,flipud(signal), ...
    signal, flipud(signal)));
ntotFrames = numel(signal);
signal_f = signal_c(ntotFrames+1:2*ntotFrames); %filtered distance from fixation cross
