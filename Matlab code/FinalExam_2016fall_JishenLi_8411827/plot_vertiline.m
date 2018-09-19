function [ output_args ] = plot_vertiline( z )
%UNTITLED3 Summary of this function goes here
%  Detailed explanation goes here
yL = get(gca,'YLim');
for i=1:length(z(2:end-1))
line([z(i+1) z(i+1)],yL,'color','g','LineStyle','--');
end 
line([z(1) z(1)],yL,'color','g','LineStyle','-.');
line([z(end) z(end)],yL,'color','g','LineStyle','-.');
set(gcf,'Color',[1,1,1]); % White background

end

