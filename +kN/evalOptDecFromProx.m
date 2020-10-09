function OptDec = evalOptDecFromProx( CostDesignDiff,CostSysFailProxDiff,numSubOpt )

numComp = size( CostDesignDiff,2 );

weight = evalCriticalWeight( CostDesignDiff,CostSysFailProxDiff );

OptDec = [];
for ww = 1:length(weight)
    weight_w = weight(ww);
    weightedSum_w = CostDesignDiff + weight_w * CostSysFailProxDiff;
    OptDec_w = zeros(numSubOpt,numComp);
    [~,OptDec_w1] = min( weightedSum_w );
    OptDec_w(1,:) = OptDec_w1;
    weightedSum_w( sub2ind( [3,numComp],OptDec_w1,1:numComp ) ) = max( weightedSum_w(:) );
    for ss = 2:numSubOpt        
        [~,OptDec1_ws] = min( weightedSum_w(:) );
        OptComp1_ws = floor( OptDec1_ws/3 ) + 1;
        OptDecVal1_ws = rem( OptDec1_ws,3 ); if ~OptDecVal1_ws; OptDecVal1_ws=3; end
        OptDec_w1( OptComp1_ws ) = OptDecVal1_ws;
        OptDec_w(ss,:) = OptDec_w1;
        weightedSum_w( OptDecVal1_ws,OptComp1_ws ) = max( weightedSum_w(:) );
    end
    OptDec = [OptDec; OptDec_w];
end
OptDec = unique(OptDec,'rows');

end



function weight = evalCriticalWeight( CostDesignDiff,CostSysFailProxDiff )

numComp = size( CostDesignDiff,2 );

weight = [];
for cc = 1:numComp
    costSysFailProxDiff_c = CostSysFailProxDiff(:,cc); costDesignDiff_c = CostDesignDiff(:,cc);
    slopeFromD1_c = -( costDesignDiff_c(2:3)-costDesignDiff_c(1) ) ./ ( costSysFailProxDiff_c(2:3) - costSysFailProxDiff_c(1) );
    if slopeFromD1_c(1) > slopeFromD1_c(2) % D1-D2 slope > D1-D3 slope
        weight = [weight; slopeFromD1_c(1)];
        slopeFromD2_c = -( costDesignDiff_c(3)-costDesignDiff_c(2) ) ./ ( costSysFailProxDiff_c(3) - costSysFailProxDiff_c(2) );
        if slopeFromD2_c > 0
            weight = [weight; slopeFromD2_c];
        end
    else
        if slopeFromD1_c(2) > 0
            weight = [weight; slopeFromD1_c(2)];
        end
    end    
end
weight = unique(weight);

end