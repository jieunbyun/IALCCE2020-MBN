function [pf,cov,numSamp] = MCS_FlowNet( Graph,linkCapa,linkProb,nodeS,nodeT,flowTarget,maxMcsCov,maxMcsSamp )

numLink = size( Graph.Edges,1 );

numSamp = 0; cov = 1; numFail = 0;
while ( cov > maxMcsCov || numSamp < 10 ) && numSamp < maxMcsSamp
    numSamp = numSamp+1;
    sampCapa_ = zeros( numLink,1 );
    for ll = 1:numLink
        sampCapa_(ll) = randsample( linkCapa,1,true,linkProb );
    end
    Graph.Edges.Weight = sampCapa_;
    flow_ = maxflow( Graph,nodeS,nodeT );
    if flow_ < flowTarget
        numFail = numFail+1;
    end
    pf = numFail / numSamp;
    cov = sqrt( (1-pf)/pf/numSamp );    
end
