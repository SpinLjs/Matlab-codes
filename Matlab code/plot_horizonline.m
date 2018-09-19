function [ output_args ] = plot_horizonline( z )
%UNTITLED3 Summary of this function goes here
%  Detailed explanation goes here
xL = get(gca,'XLim');
line(xL,[z z],'color','g','LineStyle','-');
end

