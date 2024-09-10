function screen2png(filename, h)
%SCREEN2PNG Generate a PNG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2JPEG('filename') saves the current figure to the
% PNG file "filename".
%
% SCREEN2JPEG('filename', h) saves the figure (h) to the
% PNG file "filename".
%
% 2014-4-8 DS added 2nd input
% 2020-06-04 renamed from screen2pngDS

if nargin < 1
error('Not enough input arguments!')
end

if nargin<2
    h = gcf;
end
set(h,'InvertHardcopy','off'); %31/1/22
set(0, 'currentfigure', h); %28/3/18

oldscreenunits = get(h,'Units');
oldpaperunits = get(h,'PaperUnits');
oldpaperpos = get(h,'PaperPosition');
set(h,'Units','pixels');
scrpos = get(h,'Position');
newpos = scrpos/100;
set(h,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-dpng', filename, '-r100');
drawnow
set(h,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)