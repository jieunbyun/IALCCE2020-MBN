%{
IALCCE 2020
Sioux Falls benchmark network:
Draw figure
%}

clear; close all;

fontSize_EdgeLabel = 14;
fontSize_AxisLabel = 17;
fontSize_AxisTick = 14;
fontSize_Legend = 16;
lineWidth = 1.5;
fontSize_STnode = 15;

%% 
linkPairs = load('data\SiouxFalls_net.txt');
linkPairs = sortrows(linkPairs,'ascend');
linkNo = size(linkPairs,1);
SFgraph = graph( linkPairs(:,1),linkPairs(:,2) );

nodeCoor = load( 'data\SiouxFalls_node.txt' );
nodeCoor = nodeCoor(:,2:3);
nodeCoor = nodeCoor * 0.000025; % in. -> km 
nodeNo = size( nodeCoor,1 );

nodeS = 13; nodeT = 2;

EQCoor = [(-2:2)' (2:-1:-2)'];

figure;
SFnet_fig = plot( SFgraph,'XData',nodeCoor(:,1),'YData',nodeCoor(:,2),'NodeLabel',{}, 'LineWidth',lineWidth, 'EdgeAlpha', 1,'EdgeLabel',1:linkNo );
SFnet_fig.EdgeFontSize = fontSize_EdgeLabel-2;
hold on
SFnet_fig_EQloc = plot( EQCoor(:,1),EQCoor(:,2),'ro','markerfacecolor','r','Markersize',6 );
plot( nodeCoor([nodeS nodeT],1), nodeCoor([nodeS nodeT],2), 'mp', 'markerfacecolor','m','Markersize',12)
text(nodeCoor(nodeS,1)-1, nodeCoor(nodeS,2)-.6,'\bf{Source}','FontName','times new Roman','fontsize',fontSize_STnode)
text(nodeCoor(nodeT,1)-1, nodeCoor(nodeT,2)+.9,'\bf{Terminal}','FontName','times new Roman','fontsize',fontSize_STnode)
xlabel( 'x-direction (km)','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )
ylabel( 'y-direction (km)','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )
set(gca, 'FontSize', fontSize_AxisTick,'FontName','times new roman')
legend( SFnet_fig_EQloc,'Epicenter locations','Fontsize',fontSize_Legend,'FontName','times new roman','location','southeast')
saveas(gcf,'figure/SF_graph.emf')