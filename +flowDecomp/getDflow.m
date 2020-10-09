function [D,Dflow] = getDflow( randomGraph,numNode,nodeS,nodeT,flowTarget )

randomGraphDFlow = randomGraph;
randomGraphDFlow = addedge(randomGraphDFlow,nodeT,numNode+1,flowTarget);
[D,Dflow] = maxflow( randomGraphDFlow,nodeS,numNode+1 );
dummyLinkId = find( sum(ismember(Dflow.Edges.EndNodes,numNode+1),2) );
Dflow = rmedge(Dflow,dummyLinkId);
