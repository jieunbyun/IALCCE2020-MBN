%{
IALCCE 2020
Sioux Falls benchmark network:
Connectivity 
Approximation by sampling
%}
clear;

import flowDecomp.*
import SF.*
import mbn.*
rng(1)

linkCapa = [0 1 3];
numCapa = length(linkCapa);

samplingCovTarget = .03;
%%
[linkPairs,nodeCoor,linkCoor] = summarizeSFInfo;
nodeS = 13; nodeT = 2;
EQCoor = [(-2:2)' (2:-1:-2)'];
flowTarget = 4;

numLink = size(linkPairs,1);
numNode = size( nodeCoor,1 );

%% BN quantification
numState_M = 5; numState_L = 5; numState_S = 3;
state_S = [1.4 1.0]; % Seismic capacity of roads (PGA), g(m/s2)
stateDrop_S = .2; % drop by damage status getting worse
probVec_I = [.9 .08 .02 .1 .8 .1 .02 .08 .9]';

[CPM, var, B, State] = quantCPM_SF(numLink,linkCoor,EQCoor,linkCapa,numState_M,numState_L,numState_S,state_S,stateDrop_S,probVec_I);

%% P(Xn|m5)
linkProb = evalLinkProb_CondM5( CPM,var,B );

%% Sampling
[eventMatSys,probSamp,numSamp,probSysFail,samplingCov] = sampleFlowNetBy(linkProb,linkCapa,linkPairs,nodeS,nodeT,flowTarget,samplingCovTarget);
CPM{var.System}.C = eventMatSys; CPM{var.System}.p = ones( size(eventMatSys,1),1 ); 
numSysRule = size( eventMatSys,1 );

%% Inference
observSevereID = 22; CimID1 = observSevereID; CimID2 = 31;
CPM_condObservSevere = conds(CPM,var.I(observSevereID),3,B);
CPM_absorb_IS_condObservSevere = absorb_IS_from(CPM_condObservSevere, var,B);

CPM_absorb_IS = absorb_IS_from(CPM, var,B);

eventMat_ML = [repmat( CPM{var.M}.C,numState_L,1 ) repelem( CPM{var.L}.C,numState_M,1 )];
numState_ML = numState_M*numState_L;

mapX2Sys = cell( numLink,1 );
CPM_absorb_IS_condML = conds( CPM_absorb_IS,[var.M var.L],[1 1],B );
for cc = 1:numLink
    Xind_c = CpmId( var.X(cc),CPM_absorb_IS_condML );
    mapX2Sys{cc} = mapping( CPM_absorb_IS_condML{Xind_c}.C(:,1),eventMatSys(:,cc+1),B );
end

probVecSys = zeros( numSysRule,1 );
probVecSys_observSevere = zeros( numSysRule,1 );
for ml = 1:numState_ML
    m_ = eventMat_ML(ml,1); l_ = eventMat_ML(ml,2);
    %
    CPM_absorb_IS_condML_ml = conds( CPM_absorb_IS,[var.M var.L],[m_ l_],B );
    CPM_absorb_IS_condObservSevere_condML_ml = conds( CPM_absorb_IS_condObservSevere,[var.M var.L],[m_ l_],B );
    logProb_ml = log( CPM{var.M}.p(m_) ) + log( CPM{var.L}.p(l_) );
    for rr = 1:numSysRule
        logProb_r_ml = log(1);
        logProb_r_observeSevere_ml = log(1);
       for cc = 1:numLink
           Xind_c = CpmId( var.X(cc),CPM_absorb_IS_condML_ml );
           mapX2Sys_cr = mapX2Sys{cc}(:,rr);
           %
           logProb_r_ml = logProb_r_ml + log( sum( CPM_absorb_IS_condML_ml{Xind_c}.p( mapX2Sys_cr ) ) );
           logProb_r_observeSevere_ml = logProb_r_observeSevere_ml + log( sum( CPM_absorb_IS_condObservSevere_condML_ml{Xind_c}.p( mapX2Sys_cr ) ) );             
       end
       prob_r_ml = exp( logProb_ml+logProb_r_ml ); probVecSys(rr) = probVecSys(rr) + prob_r_ml;
       prob_r_observeSevere_ml = exp( logProb_ml+logProb_r_observeSevere_ml ); probVecSys_observSevere(rr) = probVecSys_observSevere(rr) + prob_r_observeSevere_ml;
       
       if ~rem(rr,1e3)
          disp( ['Conditioning ' num2str(ml) '; Rule ' num2str(rr) ' done'] ) 
       end
    end
    
end

%%
ruleSysFail = (eventMatSys(:,1) == 2);
probVecSys_weight = probVecSys ./ probSamp;
probSysFail = 1/numSamp * sum( probVecSys_weight(ruleSysFail) );
probSysFail_var = 1/numSamp * ( sum( (probVecSys_weight(ruleSysFail)-probSysFail).^2 ) + + sum(~ruleSysFail)*probSysFail^2 );
probSysFail_cov = sqrt(probSysFail_var/ numSamp) / probSysFail ;

ruleXind1_fail = (~mapX2Sys{CimID1}(1,:) ); ruleXind1_fail = ruleXind1_fail';
probSys_X0_ID1 = sum( probVecSys_weight(ruleSysFail&ruleXind1_fail) ) / sum( probVecSys_weight(ruleSysFail) );
probSys_X0_ID1_var = sum( probVecSys_weight(ruleSysFail).^2 .* (ruleXind1_fail(ruleSysFail)-probSys_X0_ID1).^2 ) / ( sum(probVecSys_weight(ruleSysFail)) )^2;
probSys_X0_ID1_cov = sqrt( probSys_X0_ID1_var ) / probSys_X0_ID1;

probVecSys_observSevere_weight = probVecSys_observSevere ./ probSamp;
probSys_observSevere_ID1 = sum( probVecSys_observSevere_weight(ruleSysFail&ruleXind1_fail) ) / sum( probVecSys_observSevere_weight(ruleSysFail) );
probSys_observSevere_ID1_var = sum( probVecSys_observSevere_weight(ruleSysFail).^2 .* (ruleXind1_fail(ruleSysFail)-probSys_observSevere_ID1).^2 ) / ( sum(probVecSys_observSevere_weight(ruleSysFail)) )^2;
probSys_observSevere_ID1_cov = sqrt( probSys_observSevere_ID1_var ) / probSys_observSevere_ID1;

ruleXind2_fail = (~mapX2Sys{CimID2}(1,:) ); ruleXind2_fail = ruleXind2_fail';
probSys_X0_ID2 = sum( probVecSys_weight(ruleSysFail&ruleXind2_fail) ) / sum( probVecSys_weight(ruleSysFail) );
probSys_X0_ID2_var = sum( probVecSys_weight(ruleSysFail).^2 .* (ruleXind2_fail(ruleSysFail)-probSys_X0_ID2).^2 ) / ( sum(probVecSys_weight(ruleSysFail)) )^2;
probSys_X0_ID2_cov = sqrt( probSys_X0_ID2_var ) / probSys_X0_ID2;

disp( ['P(Sys=Fail) : ' num2str( probSysFail ) ' | c.o.v. : ' num2str( probSysFail_cov )] )
disp( ['P(X_ind1=Fail | Sys=Fail) : ' num2str( probSys_X0_ID1 ) ' | c.o.v. : ' num2str( probSys_X0_ID1_cov )] )
disp( ['P(X_ind2=Fail | Sys=Fail) : ' num2str( probSys_X0_ID2 ) ' | c.o.v. : ' num2str( probSys_X0_ID2_cov )] )
disp( ['P(X_ind1=Fail | Sys=Fail, I_ind1=Severe Damage) : ' num2str( probSys_observSevere_ID1 ) ' | c.o.v. : ' num2str( probSys_observSevere_ID1_cov )] )

save demoSF_sampling
