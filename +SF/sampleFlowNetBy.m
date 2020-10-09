function [eventMatSys,probSamp,numSamp,probSysFail,samplingCov] = sampleFlowNetBy(linkProb,linkCapa,linkPairs,nodeS,nodeT,flowTarget,samplingCovTarget)

rng(1)

numLink = size(linkPairs,1);
numCapa = size(linkProb,2);

graphSampling = graph( linkPairs(:,1),linkPairs(:,2) );

numSamp = 0; samplingCov = 1e2; 
probSamp = [];
numFail = 0;
eventMatSys = [];

while samplingCov > samplingCovTarget
    numSamp = numSamp+1;
    
    sampX_ = zeros( 1,numLink );
    sampXind_ = zeros( 1,numLink );
    for cc = 1:numLink
        sampXind_(cc) = randsample( 1:numCapa,1,'true',linkProb(cc,:) );
        sampX_(cc) = linkCapa( sampXind_(cc) );
    end
    probSampX_ = exp( sum( log( linkProb( sub2ind( [numLink,numCapa],1:numLink,sampXind_ ) ) ) ) );
    
    graphSF_ = graphSampling;
    graphSF_.Edges.Weight = sampX_';
    
    maxFlow_ = maxflow( graphSF_,nodeS,nodeT );
    
    if maxFlow_ < flowTarget
        numFail = numFail+1;
        eventMatSys = [eventMatSys; 2 (numCapa+1)-sampXind_];
    else
        eventMatSys = [eventMatSys; 1 (numCapa+1)-sampXind_];
    end
    
    probSamp = [probSamp; probSampX_];
    
    probSysFail = numFail / numSamp;    
    if probSysFail && numSamp > 10      
        samplingVar = (1-probSysFail)*probSysFail / numSamp;
        samplingCov = sqrt( samplingVar ) / probSysFail;
    end
    
    if ~rem( numSamp,5e2 )
        disp( ['No. of Samples: ' num2str(numSamp) ' | System Failure Prob: ' num2str(probSysFail) ' | c.o.v.: ' num2str(samplingCov)] )
    end
end
