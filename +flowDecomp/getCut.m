function cut = getCut( GraphLow,Dflow,nodeS,nodeT,linkNodePairSet )

% % DflowLinkPairs = Dflow.Edges.EndNodes;
% % DflowLinkPairs = sort(DflowLinkPairs,2);
% % DflowOnLink = Dflow.Edges.Weight;
% % 
numNode = max( linkNodePairSet(:) );
% % 
% % GraphLowLinkPairs = GraphLow.Edges.EndNodes;
% % GraphLowFlowOnLink = GraphLow.Edges.Weight;
% % for ll = 1:size(DflowLinkPairs,1)
% %     linkID_l = find( ismember( GraphLowLinkPairs,DflowLinkPairs(ll,:),'rows' ) );
% % %     if event( linkID_l,1 ) == event( linkID_l,2 )
% % %         GraphLowFlowOnLink( linkID_l ) = GraphLowFlowOnLink( linkID_l )+1;
% % %     end
% %     GraphLowFlowOnLink( linkID_l ) = GraphLowFlowOnLink( linkID_l ) - DflowOnLink(ll);
% % end

% sameUpLowLinkId = (event( :,1 ) == event( :,2 ));
% GraphLowFlowOnLink(sameUpLowLinkId) = max(GraphLowFlowOnLink)+1;

% % GraphLow = rmedge( GraphLow,GraphLowLinkPairs(GraphLowFlowOnLink<=0,1),GraphLowLinkPairs(GraphLowFlowOnLink<=0,2) );
    
% % nodeConnectedToS = neighbors( GraphLow,nodeS );
% % nodeConnectedToS = [nodeConnectedToS; nodeS];
% % nodeDifconToS = setdiff( 1:numNode,nodeConnectedToS );
% % cut = getTraversingLinkSet( nodeConnectedToS,nodeDifconToS,linkNodePairSet );

[~,~,nodeConnectedToS] = maxflow( GraphLow,nodeS,nodeT );
nodeConnectedToS = [nodeConnectedToS; nodeS];
nodeDifconToS = setdiff( 1:numNode,nodeConnectedToS );
cut = getTraversingLinkSet( nodeConnectedToS,nodeDifconToS,linkNodePairSet );

function LinkSet = getTraversingLinkSet( nodeSet1,nodeSet2,linkNodePairSet )

LinkSet = find( ismember( linkNodePairSet(:,1),nodeSet1 ) & ismember( linkNodePairSet(:,2),nodeSet2 ) );
LinkSet = [LinkSet; find( ismember( linkNodePairSet(:,2),nodeSet1 ) & ismember( linkNodePairSet(:,1),nodeSet2 ) )];

