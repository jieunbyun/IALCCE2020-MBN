function linkProb = evalLinkProb_CondM5( CPM,var,B )
import mbn.*

numLink = length( var.X );
numCapa = size( B{var.X(1)},2 );
numState_L = size( B{var.L},2);

linkProb = zeros( numLink,numCapa );

CPM_absorb_IS = CPM;
for cc = 1:numLink
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.I(cc),B );
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.S(cc),B );
end

CPM_PX_m5 = CPM_absorb_IS;
CPM_PX_m5 = conds( CPM_PX_m5,var.M,5,B );
for ll = 1:numState_L
   CPM_PX_m5_l = CPM_PX_m5;
   CPM_PX_m5_l = conds( CPM_PX_m5_l,var.L,ll,B );
   for cc = 1:numLink
       CPM_PX_m5_lc = sumProductElimVar( CPM_PX_m5_l( [3+cc var.L] ),var.L,B );
       linkProb( cc,: ) = linkProb( cc,: ) + CPM_PX_m5_lc{1}.p(end:-1:1)';
   end 
end