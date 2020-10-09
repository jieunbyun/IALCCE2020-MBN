function CPM_absorb_IS = absorb_IS_from(CPM, var,B)

import mbn.*
numLink = length( var.X );

CPM_absorb_IS = CPM;
for cc = 1:numLink
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.I(cc),B );
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.S(cc),B );
end