function [CPM, var, B, State] = quantCPM_3OutOfN_DM(numComp,compDamageProb,compDamageProbInc,Rel_Comp,RelDrop_Comp,Cost)

import mbn.*
var_ = 0;

% Dn
for cc = 1:numComp
    var_ = var_+1;
    var.D(cc,1) = var_;
    CPM{var_,1} = cpmjoint( var_,(1:3)',[1 1 1]' );
    State{var_,1} = {'Best' 'Middle' 'Cheap'}'; % choice on component type
    B{var_,1} = eye(3);
end

% Sn|Dn
eventMat_Sn = [repmat((1:3)',3,1) repelem((1:3)',3,1)];
compFailProb2 = compDamageProb+compDamageProbInc;
compFailProb3 = compDamageProb+2*compDamageProbInc;
probVec_Sn = [1-3*compDamageProb 2*compDamageProb compDamageProb 1-3*compFailProb2 2*compFailProb2 compFailProb2 1-3*compFailProb3 2*compFailProb3 compFailProb3]';
for cc = 1:numComp
    var_ = var_+1;
    var.S(cc,1) = var_;
    CPM{var_,1} = cpmcond( [var_ var.D(cc)],1,eventMat_Sn,probVec_Sn );
    State{var_,1} = {'Intact' 'Slight' 'Severe'}'; % damage state
    B{var_,1} = eye(3);
end

% Xn|Cn
EventMat_Xn = [repmat([1 2]',3,1) repelem((1:3)',2,1)];
ProbVec_Xn = [Rel_Comp 1-Rel_Comp Rel_Comp-RelDrop_Comp 1-(Rel_Comp-RelDrop_Comp) Rel_Comp-2*RelDrop_Comp 1-(Rel_Comp-2*RelDrop_Comp)]';
for cc = 1:1:numComp
    var_ = var_+1;
    var.X(cc,1) = var_;
    CPM{var_,1} = cpmcond( [var_ var.S(cc)],1,EventMat_Xn,ProbVec_Xn );
    State{var_,1} = {'Survive' 'Fail'}';
    B{var_,1} = [eye(2); 1 1];
end

% Xbar_n|Xbar_(n-1),Xn
cc = 1;
var_ = var_+1;
var.Xbar(cc,1) = var_;
CPM{var_,1} = cpmcond( [var_ var.X(cc)],1,[3 2;4 1],[1 1]' );
State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1'}'; 
B{var_,1} = eye(4);

cc = cc+1;
var_ = var_+1;
var.Xbar(cc,1) = var_;
CPM{var_,1} = cpmcond( [var_ var.Xbar(cc-1,1) var.X(cc)],1,[3 3 2; 4 3 1; 4 4 2; 5 4 1],[1 1 1 1]' );
State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1' 'Surv. Comp. No. 2'}';
B{var_,1} = eye(5);

cc = cc+1;
var_ = var_+1;
var.Xbar(cc,1) = var_;
CPM{var_,1} = cpmcond( [var_ var.Xbar(cc-1,1) var.X(cc)],1,[3 3 2; 4 3 1; 4 4 2; 5 4 1; 5 5 2; 1 5 1],[1 1 1 1 1 1]' );
State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1' 'Surv. Comp. No. 2'}';
B{var_,1} = eye(5);

EventMat_Xbar_n = [3 3 2; 4 3 1; 4 4 2; 5 4 1; 5 5 2; 1 5 1; 1 1 3];
ProbVec_Xbar_n = ones(7,1);
for cc = 4:(numComp-3)
    var_ = var_+1;
    var.Xbar(cc,1) = var_;
    CPM{var_,1} = cpmcond( [var_ var.Xbar(cc-1,1) var.X(cc)],1,EventMat_Xbar_n,ProbVec_Xbar_n );
    State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1' 'Surv. Comp. No. 2'}';
    B{var_,1} = eye(5);
end

cc = cc+1;
var_ = var_+1;
var.Xbar(cc,1) = var_;
CPM{var_,1} = cpmcond( [var_ var.Xbar(cc-1,1) var.X(cc)],1,[2 3 2; 2 3 1; 2 4 2; 3 4 1; 3 5 2; 1 5 1; 1 1 3],ones(7,1) );
State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1' 'Surv. Comp. No. 2'}';
B{var_,1} = eye(5);

cc = cc+1;
var_ = var_+1;
var.Xbar(cc,1) = var_;
CPM{var_,1} = cpmcond( [var_ var.Xbar(cc-1,1) var.X(cc)],1,[2 3 2; 3 3 1; 3 4 2; 1 4 1; 1 1 3; 2 2 3],ones(6,1) );
State{var_,1} = {'System survival' 'System failure' 'Surv. Comp. No. 0' 'Surv. Comp. No. 1'}';
B{var_,1} = eye(4);

% X_(N+1)|Xbar_(N-1),XN
var_ = var_+1;
var.System = var_;
CPM{var_,1} = cpmcond( [var_ var.Xbar(numComp-1) var.X(numComp)],1,[2 3 2; 1 3 1; 1 1 3; 2 2 3],[1 1 1 1]' );
State{var_,1} = {'Survival' 'Failure'}';
B{var_,1} = eye(2);

% Cn|Dn
EventMat_Cn = [Cost ones(3,1)];
ProbVec_Cn = ones(3,1);
for cc = 1:1:numComp
    var_ = var_+1;
    var.C(cc,1) = var_;
    CPM{var_,1} = cpmcond( [var_ var.D(cc)],1,EventMat_Cn,ProbVec_Cn );
    State{var_,1} = {'Design cost'}';
    B{var_,1} = eye(3);
end

% C_(N+1)|X_(N+1)
var_ = var_+1;
var.Csys = var_;
CPM{var_} = cpmcond( [var_ var.System],1,[0 1;1 2],[1 1]' );
State{var_,1} = {'Cost for system performance'}';
B{var_,1} = eye(2);