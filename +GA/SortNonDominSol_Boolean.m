function popSortedBoolean = SortNonDominSol_Boolean( popObjVal1,popObjVal2 )

noPop = length(popObjVal1);
yetDecide = 0; nonDomin = 1; domin = -1; underInsp = 3;

popSortedBoolean = yetDecide * ones( noPop,1 );
while sum( ~popSortedBoolean )
    popIdx_ = find( ~popSortedBoolean , 1 );
    objVal1_ = popObjVal1( popIdx_ ); objVal2_ = popObjVal2( popIdx_ );
    
    popSortedBoolean( popIdx_ ) = underInsp;
    popYetDecided_ = find( ~popSortedBoolean );
    if CheckIfDominSolExist( [objVal1_ objVal2_],[popObjVal1(popYetDecided_) popObjVal2(popYetDecided_)] )
        popSortedBoolean( popIdx_ ) = domin;
    else
        popSortedBoolean( popIdx_ ) = nonDomin;
        dominatedSolBoolean = CheckDominatedSol( [objVal1_ objVal2_],[popObjVal1(popYetDecided_) popObjVal2(popYetDecided_)] );
        popSortedBoolean( popYetDecided_(dominatedSolBoolean) ) = domin;
    end    
end

popSortedBoolean = (popSortedBoolean>0);
    
end
    
function dominSolCheck = CheckIfDominSolExist( popObjVals,otherObjVals )
    noObj = length(popObjVals);
    noOthers = size( otherObjVals,1 );
    
    dominCheckVector = ones( noOthers,1 );
    for objVal = 1:noObj
        dominCheckVector = dominCheckVector .* ( otherObjVals(:,objVal)<popObjVals(objVal) );
    end
    
    dominSolCheck = ( sum( dominCheckVector ) > 0);
end

function dominatedSolBoolean = CheckDominatedSol( popObjVals,otherObjVals )
    noObj = length(popObjVals);
    noOthers = size( otherObjVals,1 );
    
    dominCheckVector = ones( noOthers,1 );
    for objVal = 1:noObj
        dominCheckVector = dominCheckVector .* ( otherObjVals(:,objVal)>popObjVals(objVal) );
    end
    
    dominatedSolBoolean = ( dominCheckVector > 0);
end