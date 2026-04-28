function [ figformat ] = savefigure( name )
% name: Name of the figure, example: 'figure1'
% tupe: Format of the figure, example: 'fig','espc'
% If no input is provided, figure saved twice with fig and espc

ni=nargin;
set(gcf,'PaperPositionMode','auto') 
print -depsc -r0 pracfigsave.eps 
if ni==1
    saveas(gcf,name,'fig') 
    saveas(gcf,name,'epsc')
%     saveas(gcf,'Results/BIIL_8_3_2011/Figures/cf_lag.eps','epsc')
    figformat{1}='fig';figformat{2}='epsc';
%     keyboard
else
    saveas(gcf,name,type) 
    figformat{1}=type;
end

end

