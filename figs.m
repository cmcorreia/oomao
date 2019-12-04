function out = figs(h, filename)
try cd(['_fig_' filename])
catch
    mkdir(['_fig_' filename])
    cd(['_fig_' filename])
end
savefig([filename]);
print(filename,'-depsc','-tiff', '-painters', '-r864');
print(filename,'-dpng');

fprintf('Printing: %s to %s\n',filename, pwd)
cd ..
end
