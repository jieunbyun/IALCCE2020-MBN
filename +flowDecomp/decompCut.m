function [eventFail,eventFailProb,eventUnspecSetIn,eventUnspecSetInProb] = decompCut( event,cut,cutCapaId,linkProb )

import flowDecomp.*

bndLow = event(:,1); bndUp = event(:,2);

eventUnspecSetIn = {};
for ll = 1:length(cut)
    cut_l = cut(ll);
    cutCapaId_l = cutCapaId(ll);
    bndLowOut = bndLow;
    bndLowOut(cut_l) = cutCapaId_l+1;
    eventOut = [bndLowOut bndUp];
    eventUnspecSetIn = [eventUnspecSetIn; {eventOut}];
    bndUp(cut_l) = cutCapaId_l;
end
eventFail = [bndLow bndUp];

eventFailProb = getEventSetProb( {eventFail},linkProb );
eventUnspecSetInProb = getEventSetProb( eventUnspecSetIn,linkProb );