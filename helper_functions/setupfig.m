function [] = setupfig(wid,ht,ft)
set(gcf ,'windowStyle','normal');
set(gcf,'color','w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 wid ht]);

set(findall(gcf,'-property','FontSize'),'FontSize',ft)