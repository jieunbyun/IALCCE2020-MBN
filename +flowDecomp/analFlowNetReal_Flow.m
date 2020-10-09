function [eventSurvSet, eventSurvSetProb, eventFailSet, eventFailSetProb, eventUnspecSet, eventUnspecSetProb, numFlowAnal] = ...
    analFlowNetReal_Flow( Graph,linkCapa,linkProb,nodeS,nodeT,flowTarget,numNode,bndWidth,maxNumSet )

import flowDecomp.*

numLink = size( Graph.Edges,1 );
linkNodePairSet = Graph.Edges.EndNodes;

eventUnspecSet = {[ones(numLink,1) length(linkCapa)*ones(numLink,1)]}; % bnd: [lower, upper]
eventUnspecSetProb = 1;

eventSurvSet = {}; eventSurvSetProb = [];
eventFailSet = {}; eventFailSetProb = [];
numFlowAnal = 0;
countForMessage = 1;
while sum(eventUnspecSetProb)/sum(eventFailSetProb) > bndWidth && length( [eventSurvSetProb;eventFailSetProb;eventUnspecSetProb] ) < maxNumSet
    [~,eventId] = max(eventUnspecSetProb);
    event = eventUnspecSet{eventId};
    eventProb = eventUnspecSetProb(eventId);
    
    Graph.Edges.Weight = linkCapa(event(:,2))';
    [D,Dflow] = getDflow( Graph,numNode,nodeS,nodeT,flowTarget );
    numFlowAnal = numFlowAnal+1;

    if D < flowTarget
        eventFailSet = [eventFailSet; {event}];
        eventFailSetProb = [eventFailSetProb; eventProb];
        eventUnspecSetIn = {}; eventUnspecSetInProb = [];
    else
        Graph.Edges.Weight = linkCapa(event(:,1))';
        flowLow = maxflow( Graph,nodeS,nodeT );
        numFlowAnal = numFlowAnal+1;
        
        if flowLow < flowTarget
        
            [eventSurv,eventSurvProb,eventUnspecSetIn,eventUnspecSetInProb] = decompFlow( event,Dflow,linkNodePairSet,linkCapa,linkProb );
            
            eventSurvSet = [eventSurvSet; {eventSurv}];
            eventSurvSetProb = [eventSurvSetProb; eventSurvProb];
            
        else
            eventSurvSet = [eventSurvSet; {event}];
            eventSurvSetProb = [eventSurvSetProb; eventProb];
            eventUnspecSetIn = {}; eventUnspecSetInProb = [];
        end
    end
    
    [eventUnspecSet,eventUnspecSetProb] = updateUnspecSet( eventUnspecSet,eventUnspecSetProb,eventId,eventUnspecSetIn,eventUnspecSetInProb );
    
    numSetAll = length(eventUnspecSetProb)+length(eventFailSetProb)+length(eventSurvSetProb);
    if numSetAll > 5e2*countForMessage
        countForMessage = countForMessage+1;
        disp( ['No. of set: ' num2str(numSetAll) ' | Unspec. Prob / Fail. Prob: ' num2str( sum(eventUnspecSetProb)/sum(eventFailSetProb) )] )
    end
    
end