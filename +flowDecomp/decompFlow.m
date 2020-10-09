function [eventSurv,eventSurvProb,eventUnspecSetIn,eventUnspecSetInProb] = decompFlow( event,Dflow,linkNodePairSet,linkCapa,linkProb )

import flowDecomp.*

bndLow = event(:,1); bndUp = event(:,2);
DflowWeight = Dflow.Edges.Weight;
DflowLinkPair = Dflow.Edges.EndNodes;

linkID = []; flow_diff = []; capaLowId = [];
for ll = 1:length(DflowWeight)
    linkpair_l = DflowLinkPair(ll,:);
    linkpair_l = sort(linkpair_l);
    linkflow_l = DflowWeight(ll);
    linkId_l = find( ismember(linkNodePairSet,linkpair_l,'rows') );
    
    minFlow_l = linkCapa(event(linkId_l,1 ));
    if minFlow_l < linkflow_l
        linkID = [linkID; linkId_l];
        capaLowId_l = find( linkCapa < linkflow_l , 1, 'last' );
        capaLowId = [capaLowId; capaLowId_l];
        flow_diff = [flow_diff; linkCapa(capaLowId_l)-minFlow_l];
    end
end

[~,flowSort] = sort( flow_diff,'descend' );
linkID = linkID(flowSort); capaLowId = capaLowId(flowSort);

eventUnspecSetIn = {};
for ll = 1:length(linkID)
    linkId_l = linkID(ll);
    capaLowId_l = capaLowId(ll);

    bndUpOut = bndUp;
    bndUpOut(linkId_l) = capaLowId_l;
    eventOut = [bndLow bndUpOut];
    eventUnspecSetIn = [eventUnspecSetIn; {eventOut}];
    bndLow(linkId_l) = capaLowId_l+1;

end
eventSurv = [bndLow bndUp];

eventSurvProb = getEventSetProb( {eventSurv},linkProb );
eventUnspecSetInProb = getEventSetProb( eventUnspecSetIn,linkProb );
