%{
IALCCE 2020
3-out-of-N system
Exact inference
%}
clear;

import mbn.*
import kN.*
%% Quantification
numComp = 30;
probCond = [.8 .15 .05]'; % Cn: [Intact SlightDamage SevereDamage]
Rel_Comp = .5; RelDrop_Comp = .2; % Xn: P(Xn=Surv | Cn = Intact)=Rel_Comp; P(Xn=Surv | Cn = SlightDamage) = Rel_Comp-RelDrop_Comp --> P(Xn=Surv | Cn = SeverDamage) = Rel_Comp-2*RelDrop_Comp
InspError = .05; % In: P(In==Cn+-1)=2*InspError; P(In==Cn+-2) = InspError
[CPM, var, B, State] = quantCPM_3OutOfN(numComp,probCond,Rel_Comp,RelDrop_Comp,InspError);

%% Inference
% [1] System reliability
varElimOrder = [var.I; var.S];
for cc = 1:(numComp-1)
    varElimOrder = [varElimOrder; var.X(cc); var.Xbar(cc)]; 
end
varElimOrder = [varElimOrder; var.X(numComp)];

CPM_sysRel = CPM;
for vv = 1:length(varElimOrder)
    varElim_v = varElimOrder(vv);
    CPM_sysRel = sumProductElimVar( CPM_sysRel,varElim_v,B );
end

disp( ['P(System Failure) = ' num2str( CPM_sysRel{1}.p(1) )] )

%% Inference by Evidence
Observ = [ones(numComp/3,1); 2*ones(numComp/3,1); 3*ones(numComp/3,1)];
CPM_observ = conds( CPM,var.I,Observ,B );
for vv = 1:length(varElimOrder)
    varElim_v = varElimOrder(vv);
    CPM_observ = sumProductElimVar( CPM_observ,varElim_v,B );
end
CPM_observ{1}.p = exp( log(CPM_observ{1}.p) - log( sum(CPM_observ{1}.p) ) );
disp( ['P(System Failure | Observ. ) = ' num2str( CPM_observ{1}.p(1) )] )


%% Component Importance Measure P(XN=fail | Sys=fail)
CPM_CIM_observ1 = conds( CPM,[var.I; var.System],[Observ; 2],B );
varElimOrder_CIM_observ1 = setdiff(varElimOrder,var.X(1),'stable');
varElimOrder_CIM_observ1 = [varElimOrder_CIM_observ1; var.System];
for vv = 1:length(varElimOrder_CIM_observ1)
    varElim_v = varElimOrder_CIM_observ1(vv);
    CPM_CIM_observ1 = sumProductElimVar( CPM_CIM_observ1,varElim_v,B );
end
CPM_CIM_observ1{2}.p = exp( log(CPM_CIM_observ1{2}.p)-log(sum(CPM_CIM_observ1{2}.p)) );

CPM_CIM_observ2 = conds( CPM,[var.I; var.System],[Observ; 2],B );
varElimOrder_CIM_observ2 = setdiff(varElimOrder,var.X(11),'stable');
varElimOrder_CIM_observ2 = [varElimOrder_CIM_observ2; var.System];
for vv = 1:length(varElimOrder_CIM_observ2)
    varElim_v = varElimOrder_CIM_observ2(vv);
    CPM_CIM_observ2 = sumProductElimVar( CPM_CIM_observ2,varElim_v,B );
end
CPM_CIM_observ2{2}.p = exp( log(CPM_CIM_observ2{2}.p)-log(sum(CPM_CIM_observ2{2}.p)) );

CPM_CIM_observ3 = conds( CPM,[var.I; var.System],[Observ; 2],B );
varElimOrder_CIM_observ3 = setdiff(varElimOrder,var.X(21),'stable');
varElimOrder_CIM_observ3 = [varElimOrder_CIM_observ3; var.System];
for vv = 1:length(varElimOrder_CIM_observ3)
    varElim_v = varElimOrder_CIM_observ3(vv);
    CPM_CIM_observ3 = sumProductElimVar( CPM_CIM_observ3,varElim_v,B );
end
CPM_CIM_observ3{2}.p = exp( log(CPM_CIM_observ3{2}.p)-log(sum(CPM_CIM_observ3{2}.p)) );

disp( ['P(Xn Fail | In=Intact, System Failure ) = ' num2str( CPM_CIM_observ1{2}.p(1) )] )
disp( ['P(Xn Fail | In=Slight Damage, System Failure ) = ' num2str( CPM_CIM_observ2{2}.p(1) )] )
disp( ['P(Xn Fail | In=Severe Damage, System Failure ) = ' num2str( CPM_CIM_observ3{2}.p(1) )] )