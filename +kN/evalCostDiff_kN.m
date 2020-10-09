function [CostDesignDiff,CostSysFailProxDiff] = evalCostDiff_kN( baseD,Cost,CPM_noCost,var,B,varElimOrder )

import kN.*
import mbn.*

numComp = length( baseD );

[CostDesign_baseD,CostSysFail_baseD] = evalPopKn( baseD,Cost,CPM_noCost,var,B,varElimOrder );
CostDesign = zeros( 3,numComp ); CostSysFailProx = zeros( 3,numComp );
for cc = 1:numComp
    varElimOrder_c = setdiff( varElimOrder,[var.D(cc) var.S(cc) var.X(cc)] );
    
    CPM_noCost_c = CPM_noCost;
    CPM_noCost_c = conds( CPM_noCost_c,var.D( setdiff(1:numComp,cc) ),baseD( setdiff(1:numComp,cc) ),B );
    for var_ = varElimOrder_c'
        CPM_noCost_c = sumProductElimVar( CPM_noCost_c,var_,B );
    end    
    
    for dd = 1:3
        CPM_noCost_c_d = CPM_noCost_c;
        CPM_noCost_c_d = conds( CPM_noCost_c_d,var.D(cc),dd,B );
        varElimOrder_d = [var.D(cc) var.S(cc) var.X(cc)];
        for var_ = varElimOrder_d
            CPM_noCost_c_d = sumProductElimVar(CPM_noCost_c_d,var_,B);
        end
        CostSysFailProx(dd,cc) = CPM_noCost_c_d{1}.p(2);
        
        baseD_d = baseD;
        baseD_d(cc) = dd;
        CostDesign(dd,cc) = sum( Cost( baseD_d ) );
    end
    
end
CostDesignDiff = CostDesign-CostDesign_baseD;
CostSysFailProxDiff = CostSysFailProx - CostSysFail_baseD;