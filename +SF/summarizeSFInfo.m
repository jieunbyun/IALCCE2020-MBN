function [linkPairs,nodeCoor,linkCoor] = summarizeSFInfo

linkPairs = load('data\SiouxFalls_net.txt');
linkPairs = sortrows(linkPairs,'ascend');

nodeCoor = load( 'data\SiouxFalls_node.txt' );
nodeCoor = nodeCoor(:,2:3);
nodeCoor = nodeCoor * 0.000025; % in. -> km 

linkCoor = .5*( nodeCoor( linkPairs(:,1),: ) + nodeCoor( linkPairs(:,2),: ) );