%{
IALCCE 2020
Sioux Falls benchmark network:
Connectivity 
Deterministic approximation
%}

import flowDecomp.*
import SF.*
import mbn.*

linkCapa = [0 1 3];
numCapa = length(linkCapa);

bndWidth = 0.05; maxNumSet = 1e5;
%% 
linkPairs = load('data\SiouxFalls_net.txt');
linkPairs = sortrows(linkPairs,'ascend');
numLink = size(linkPairs,1);
SFgraph = graph( linkPairs(:,1),linkPairs(:,2) );

nodeCoor = load( 'data\SiouxFalls_node.txt' );
nodeCoor = nodeCoor(:,2:3);
nodeCoor = nodeCoor * 0.000025; % in. -> km 
numNode = size( nodeCoor,1 );

linkCoor = .5*( nodeCoor( linkPairs(:,1),: ) + nodeCoor( linkPairs(:,2),: ) );

nodeS = 13;
nodeT = 2;

EQCoor = [(-2:2)' (2:-1:-2)'];

SF.nodeNo = numNode;
SF.linkNo = numLink;
SF.linkPairs = linkPairs;
%% BN quantification
numState_M = 5; numState_L = 5; numState_S = 3;
state_S = [1.4 1.0]; % Seismic capacity of roads (PGA), g(m/s2)
stateDrop_S = .2; % drop by damage status getting worse
probVec_I = [.9 .08 .02 .1 .8 .1 .02 .08 .9]';
[CPM, var, B, State] = quantCPM_SF(numLink,linkCoor,EQCoor,linkCapa,numState_M,numState_L,numState_S,state_S,stateDrop_S,probVec_I);

%% P(Xn|m5)
linkProb = evalLinkProb_CondM5( CPM,var,B );

%% Decomposition -- exact
Graph = graph( linkPairs(:,1),linkPairs(:,2) );

flowTarget = 4;
[DFC.eventSurvSet, DFC.eventSurvSetProb, DFC.eventFailSet, DFC.eventFailSetProb, DFC.eventUnspecSet, DFC.eventUnspecSetProb, DFC.numFlowAnal, DFC.numDecompFlow, DFC.numDecompCut] = ...
                analFlowNetReal_FlowCut( Graph,linkCapa,linkProb,nodeS,nodeT,flowTarget,numNode,bndWidth,maxNumSet );

[eventMat_System,probVec_System] = quantCPMSys_netFlow( DFC.eventSurvSet,DFC.eventFailSet );
CPM{var.System}.C = eventMat_System; CPM{var.System}.p = probVec_System; 

%% Inference
CPM_absorb_IS = CPM;
for cc = 1:numLink
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.I(cc),B );
    CPM_absorb_IS = sumProductElimVar( CPM_absorb_IS,var.S(cc),B );
end

observSevereID = 22; CimID1 = observSevereID;
CimID2 = 31;
CPM_absorb_IS_condObservSevere = CPM;
CPM_absorb_IS_condObservSevere = conds(CPM_absorb_IS_condObservSevere,var.I(observSevereID),3,B);
for cc = 1:numLink
    CPM_absorb_IS_condObservSevere = sumProductElimVar( CPM_absorb_IS_condObservSevere,var.I(cc),B );
    CPM_absorb_IS_condObservSevere = sumProductElimVar( CPM_absorb_IS_condObservSevere,var.S(cc),B );
end

eventMat_ML = [repmat( CPM{var.M}.C,numState_L,1 ) repelem( CPM{var.L}.C,numState_M,1 )];
numState_ML = numState_M*numState_L;
probSysFail = zeros(2,1);
probXFail_SysFail1 = zeros( 6,1 ); eventXFail_SysFail = [repmat( (1:numCapa)',2,1 ) repelem( (1:2)',numCapa,1 )];
probXFail_SysFail_condISevere = zeros( 6,1 ); 
probXFail_SysFail2 = zeros( 6,1 ); 
for ml = 1:numState_ML
    m_ = eventMat_ML(ml,1); l_ = eventMat_ML(ml,2);
    %
    CPM_absorb_IS_ml = conds( CPM_absorb_IS,[var.M var.L],[m_ l_],B );
    for cc = setdiff( 1:numLink,[CimID1 CimID2] )
        CPM_absorb_IS_ml = sumProductElimVar( CPM_absorb_IS_ml,var.X(cc),B );
    end
    CPM_absorb_IS_ml = sumProductElimVar( CPM_absorb_IS_ml,var.M,B );
    CPM_absorb_IS_ml = sumProductElimVar( CPM_absorb_IS_ml,var.L,B );
    
    CPM_absorb_IS_XFail_SysFail1_ml = sumProductElimVar( CPM_absorb_IS_ml,var.X(CimID2),B );
    CPM_absorb_IS_XFail_SysFail2_ml = sumProductElimVar( CPM_absorb_IS_ml,var.X(CimID1),B );
    CPM_absorb_IS_System_ml = sumProductElimVar( CPM_absorb_IS_XFail_SysFail2_ml,var.X(CimID2),B );
    
    probSysFail = probSysFail + CPM_absorb_IS_System_ml{1}.p;
    for rr = 1:size( CPM_absorb_IS_XFail_SysFail1_ml{1}.C,1 )
        c_r = CPM_absorb_IS_XFail_SysFail1_ml{1}.C(rr,:);
        [~,c_r_idx] = ismember( eventXFail_SysFail,c_r,'rows' ); c_r_idx = logical(c_r_idx);
        probXFail_SysFail1(c_r_idx) = probXFail_SysFail1(c_r_idx) + CPM_absorb_IS_XFail_SysFail1_ml{1}.p(rr);
    end
    for rr = 1:size( CPM_absorb_IS_XFail_SysFail2_ml{1}.C,1 )
        c_r = CPM_absorb_IS_XFail_SysFail2_ml{1}.C(rr,:);
        [~,c_r_idx] = ismember( eventXFail_SysFail,c_r,'rows' ); c_r_idx = logical(c_r_idx);
        probXFail_SysFail2(c_r_idx) = probXFail_SysFail2(c_r_idx) + CPM_absorb_IS_XFail_SysFail2_ml{1}.p(rr);
    end
    %
    CPM_absorb_IS_condObservSevere_ml = CPM_absorb_IS_condObservSevere;
    CPM_absorb_IS_condObservSevere_ml = conds( CPM_absorb_IS_condObservSevere_ml,[var.M var.L],[m_ l_],B );
    for cc = setdiff( 1:numLink,observSevereID )
        CPM_absorb_IS_condObservSevere_ml = sumProductElimVar( CPM_absorb_IS_condObservSevere_ml,var.X(cc),B );
    end
    CPM_absorb_IS_condObservSevere_ml = sumProductElimVar( CPM_absorb_IS_condObservSevere_ml,var.M,B );
    CPM_absorb_IS_condObservSevere_ml = sumProductElimVar( CPM_absorb_IS_condObservSevere_ml,var.L,B );
    for rr = 1:size( CPM_absorb_IS_condObservSevere_ml{1}.C,1 )
        c_r = CPM_absorb_IS_condObservSevere_ml{1}.C(rr,:);
        [~,c_r_idx] = ismember( eventXFail_SysFail,c_r,'rows' ); c_r_idx = logical(c_r_idx);
        probXFail_SysFail_condISevere(c_r_idx) = probXFail_SysFail_condISevere(c_r_idx) + CPM_absorb_IS_condObservSevere_ml{1}.p(rr);
    end
    %
    
    disp( [num2str(ml) '-th conditioning done'] )
end    

probSysFailBnd = [probSysFail(2) 1-probSysFail(1)];
probXFail_SysFail1_Bnd = [sum(probXFail_SysFail1(5:6)) 1-sum(probXFail_SysFail1(1:4))];
probXFail_SysFail2_Bnd = [sum(probXFail_SysFail2(5:6)) 1-sum(probXFail_SysFail2(1:4))];
disp( ['P(Sys=Fail) : [' num2str( probSysFailBnd(1) ) ' ' num2str( probSysFailBnd(2) ) ']'] )
disp( ['P(X_ind1=Fail | Sys=Fail) : [' num2str( probXFail_SysFail1_Bnd(1)/probSysFailBnd(2) ) ' ' num2str( probXFail_SysFail1_Bnd(2)/probSysFailBnd(1) ) ']'] )
disp( ['P(X_ind2=Fail | Sys=Fail) : [' num2str( probXFail_SysFail2_Bnd(1)/probSysFailBnd(2) ) ' ' num2str( probXFail_SysFail2_Bnd(2)/probSysFailBnd(1) ) ']'] )

CPM_IS_observ = sumProductElimVar( CPM([var.I(observSevereID), var.S(observSevereID)]),var.S(observSevereID),B );
probXFail_SysFail_condI_Bnd = [sum(probXFail_SysFail_condISevere(5:6)) 1-sum(CPM_IS_observ{1}.p(1:2))-sum(probXFail_SysFail_condISevere(1:4))];
probSysFail_condI_Bnd = [sum( probXFail_SysFail_condISevere(4:6) ) 1-sum(CPM_IS_observ{1}.p(1:2))-sum(probXFail_SysFail_condISevere(1:3))];
disp( ['P(X_ind1=Fail | Sys=Fail, I_ind1=Severe Damage) : [' num2str( probXFail_SysFail_condI_Bnd(1)/probSysFail_condI_Bnd(2) ) ' ' num2str( probXFail_SysFail_condI_Bnd(2)/probSysFail_condI_Bnd(1) ) ']'] )

save demoSF_deterministic