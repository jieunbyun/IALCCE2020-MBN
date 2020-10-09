function [cut,cutCapaId] = getCutCapaId( cut,event,linkCapa,linkProb,flowTarget )

[cutCapaId,cutProb,flowThroughCut] = initCutCapaEval( cut,event,linkCapa,linkProb,flowTarget );

numCut = length(cut);
while flowThroughCut < flowTarget && sum(cutProb)
    [~,linkIdInc] = max(cutProb);
    capaId = cutCapaId(linkIdInc);
    flowAdd = linkCapa(capaId+1) - linkCapa(capaId);
    flowThroughCut = flowThroughCut+flowAdd;
    cutCapaId(linkIdInc) = cutCapaId(linkIdInc)+1;
    for ll = 1:numCut
        if cutProb(ll)
            link_l = cut( ll );
            linkCapaId_l = cutCapaId(ll);
            if linkCapaId_l < event(link_l,2) && ...
                flowThroughCut + linkCapa( linkCapaId_l+1 )-linkCapa( linkCapaId_l ) < flowTarget
                cutProb(ll) = linkProb( link_l,linkCapaId_l+1 );
            else
                cutProb(ll) = 0;
            end
        end
    end            
end
cutIdKeep = ( cutCapaId < event(cut,2) );
cutCapaId = cutCapaId( cutIdKeep );
cut = cut( cutIdKeep );

flow_diff = linkCapa( event( cut,2 ) ) - linkCapa( cutCapaId );
[~,cutSort] = sort( flow_diff,'descend' );
cut = cut(cutSort); cutCapaId = cutCapaId(cutSort);


function [cutCapaIdInit,cutProbInit,flowThroughCut] = initCutCapaEval( cut,event,linkCapa,linkProb,flowTarget )

numCut = length(cut);
cutCapaIdInit = event( cut,1 );
cutProbInit = zeros( numCut,1 );
flowThroughCut = sum( linkCapa(cutCapaIdInit) );
for ll = 1:numCut
    link_l = cut(ll);
    linkCapaId_l = cutCapaIdInit(ll);
%     cutCapaIdInit(ll) = linkCapaId_l;
%     flowThroughCut = flowThroughCut+linkCapa( linkCapaId_l );
    if isIncPossible( linkCapaId_l,event(link_l,2),flowThroughCut,flowTarget,linkCapa )
        cutProbInit(ll) = linkProb( link_l,linkCapaId_l+1 );
    end
end

function inc = isIncPossible( linkCapaId,CapaIdUp,flowCurrent,flowMax,linkCapa )

inc = linkCapaId < CapaIdUp && ...
            flowCurrent + linkCapa( linkCapaId+1 )-linkCapa( linkCapaId ) < flowMax;