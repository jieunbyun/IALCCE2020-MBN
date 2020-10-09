function [CostDesign,CostSysFail] = evalPopKn( pop,Cost,CPM_noCost,var,B,varElimOrder )

import mbn.*

numPop = size(pop,1);

CostDesign = zeros(numPop,1); CostSysFail = zeros(numPop,1);
for pp = 1:numPop
    pop_p = pop(pp,:);
    CostDesign(pp) = sum( Cost(pop_p) );
    
    CPM_p = CPM_noCost;
    CPM_p = conds( CPM_p,var.D,pop_p,B );
    for vv = 1:length(varElimOrder)
        varElim_v = varElimOrder(vv);
        CPM_p = sumProductElimVar( CPM_p,varElim_v,B );
    end
    CostSysFail(pp) = CPM_p{1}.p(2);
end