%{
% ------------HEADER-----------------
% Purpose:
%   Save Matlab plots in *.png, *.fig and *.eps
%INPUT:
%   h            ::  handle to a figure (default gca, i.e. get current axis)
%   filename     ::  the filename without extension (see example below)
%OUTPUT:
%   null         ::  The figures in a unique folder starting with _fig_
%
%Created by      ::  Carlos Correia
%Creation date   ::  2017
%Change Record:  ::
% ------------HEADER END-----------------

Example #1

figure, plot(1:10) % a new figure with a straight line
figs(gca,'myNewFigureFilename') %

Example #2
h = figure; % create a handle 'h' to a new figure
plot(1:10) % a new figure with a straight line
figs(h,'myNewFigureFilename') %

%}

function figs(h, filename)
try cd(['_fig_' filename])
catch
    mkdir(['_fig_' filename])
    cd(['_fig_' filename])
end
savefig(filename);
print(filename,'-depsc','-tiff', '-painters', '-r864');
print(filename,'-dpng');

fprintf('Printing: %s to %s\n',filename, pwd)
cd ..
end
