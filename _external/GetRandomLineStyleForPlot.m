function [styleString] = GetRandomLineStyleForPlot()
% This function creates the random style for your plot
% Colors iterate over all colors except for white one
  markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
  lineStyles = {'-', '--', ':', '-.'};
  colors = {'y', 'm', 'c', 'r', 'g', 'b', 'k'};
  styleString = strcat(markers(randi(length(markers), 1) ), ...
                lineStyles(randi(length(lineStyles), 1) ), ...
                colors(randi(length(colors), 1) ) );

end