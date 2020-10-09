function [eventSurvSet, eventSurvSetProb, eventFailSet, eventFailSetProb, eventUnspecSet, eventUnspecSetProb, numFlowAnal] = ...
    analFlowNetReal_Cut( Graph,linkCapa,linkProb,nodeS,nodeT,flowTarget,numNode,bndWidth,maxNumSet )

import flowDecomp.*

numLink = size( Graph.Edges,1 );
linkNodePairSet = Graph.Edges.EndNodes;

eventUnspecSet = {[ones(numLink,1) length(linkCapa)*ones(numLink,1)]}; % bnd: [lower, upper]
eventUnspecSetProb = 1;

eventSurvSet = {}; eventSurvSetProb = [];
eventFailSet = {}; eventFailSetProb = [];

numFlowAnal = 0;
while sum(eventUnspecSetProb) > bndWidth && length( [eventSurvSetProb;eventFailSetProb;eventUnspecSetProb] ) < maxNumSet
    [~,eventId] = max(eventUnspecSetProb);
    event = eventUnspecSet{eventId};
    eventProb = eventUnspecSetProb(eventId);
    
    Graph.Edges.Weight = linkCapa(event(:,2))';
    D = getDflow( Graph,numNode,nodeS,nodeT,flowTarget );
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
            cut = getCut( Graph,nodeS,nodeT,linkNodePairSet );
            [cut,cutCapaId] = getCutCapaId( cut,event,linkCapa,linkProb,flowTarget );
            [eventFail,eventFailProb,eventUnspecSetIn,eventUnspecSetInProb] = decompCut( event,cut,cutCapaId,linkProb );
            
            eventFailSet = [eventFailSet; {eventFail}];
            eventFailSetProb = [eventFailSetProb; eventFailProb];
            
        else
            eventSurvSet = [eventSurvSet; {event}];
            eventSurvSetProb = [eventSurvSetProb; eventProb];
            eventUnspecSetIn = {}; eventUnspecSetInProb = [];
        end    
    end
    [eventUnspecSet,eventUnspecSetProb] = updateUnspecSet( eventUnspecSet,eventUnspecSetProb,eventId,eventUnspecSetIn,eventUnspecSetInProb );
end