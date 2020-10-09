%{
IALCCE 2020
Eastern Massachusettes network:
Draw figure
%}

clear; close all;

fontSize_AxisLabel = 17;
fontSize_AxisTick = 14;
lineWidth = 1;

load data/ema
numLink = size( linkPairs,1 );
graphDraw = graph( linkPairs(:,1),linkPairs(:,2) );

figure;
graph_fig = plot( graphDraw,'XData',nodeCoor(:,1),'YData',nodeCoor(:,2),'NodeLabel',{}, 'LineWidth',lineWidth, 'EdgeAlpha', 1 );

set(gca, 'FontSize', fontSize_AxisTick,'FontName','times new roman')
xlabel( 'x-direction (km)','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )
ylabel( 'y-direction (km)','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )

saveas(gcf,'figure/EMA_graph.emf')