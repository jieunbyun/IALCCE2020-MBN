%{
IALCCE 2020
3-out-of-N system
Optimization
%}
clear; close all;

import mbn.*
import kN.*
%% Quantification
numComp = 30;
compDamageProb = .2; % P( Sn=Slight Damage|Dn=best ) = 2*compDamageProb; P( Sn=Severe Damage|Dn=best ) = compDamageProb; 
compDamageProbInc = .05; % Dn=middle -->  compDamageProb = compDamageProb + compFailProbInc; Dn = Cheap --> compDamageProb + 2*compFailProbInc
Rel_Comp = .5; RelDrop_Comp = .2; % Xn: P(Xn=Surv | Cn = Intact)=Rel_Comp; P(Xn=Surv | Cn = SlightDamage) = Rel_Comp-RelDrop_Comp --> P(Xn=Surv | Cn = SeverDamage) = Rel_Comp-2*RelDrop_Comp
Cost = [60 30 20]';
[CPM, var, B, State] = quantCPM_3OutOfN_DM(numComp,compDamageProb,compDamageProbInc,Rel_Comp,RelDrop_Comp,Cost);

numDecisions = 3*ones(numComp,1);
varElimOrder = [var.D; var.S; var.X; var.Xbar];
CPM_noCost = CPM;
CPM_noCost([var.C; var.Csys]) = [];

import GA.*
%% Genetic algorithm
% %{
rng(1)

numPop = 50; % # of populations at each gen.
numGen = 1500; % # of gen's to be examined
ratMut = .4; % Ratio of mutated populations
ratCross = .4; % Ratio of pupoluations being cross-over
probGeneChange = .15; % Ratio of mutation (genetic changes)

numPopMut = floor( numPop*ratMut ); numPopCross = floor( numPop*ratCross ); numPopRand = numPop-numPopMut-numPopCross;
numSol = zeros(numGen,1);

% First Gen.
pop = GenRand( numPop,numDecisions );
[CostDesign,CostSysFail] = evalPopKn( pop,Cost,CPM_noCost,var,B,varElimOrder );
% [pop,CostDesign,CostSysFail] = SortNonDominSol( pop,CostDesign,CostSysFail );
popSortedBoolean = SortNonDominSol_Boolean( CostDesign,CostSysFail );
pop = pop(popSortedBoolean,:); CostDesign = CostDesign(popSortedBoolean); CostSysFail = CostSysFail(popSortedBoolean);
numSol(1) = size(pop,1);

% Iteration
for gg = 2:numGen
    popNew = GenRand(numPopRand,numDecisions);
    popNew = [popNew; GenMut( numPopMut,numDecisions,pop,probGeneChange ) ];
    popNew = [popNew; GenCro( numPopCross,pop )];
    
    [CostDesignNew,CostSysFailNew] = evalPopKn( popNew,Cost,CPM_noCost,var,B,varElimOrder );
    
    pop = [pop; popNew]; CostDesign = [CostDesign; CostDesignNew]; CostSysFail = [CostSysFail; CostSysFailNew];
    pop = sort(pop,2); % symmetry between component events
    [pop,popUniqueIdx] = unique(pop,'rows'); CostDesign = CostDesign(popUniqueIdx); CostSysFail = CostSysFail(popUniqueIdx);
    
    popSortedBoolean = SortNonDominSol_Boolean( CostDesign,CostSysFail );
    pop = pop(popSortedBoolean,:); CostDesign = CostDesign(popSortedBoolean); CostSysFail = CostSysFail(popSortedBoolean);
    
    numSol(gg) = size(pop,1);
    
    if ~rem(gg,5e2)
        figure;
        plot( CostSysFail,CostDesign,'sq' ); grid on; drawnow;
        save( ['demoKN_optim_gen' num2str(gg)] )
    end
end
%}

%% Byun and Song (Under review)
% %{
import CD.*

rng(1)
numIter = 20;
numSubOpt = 5;

BaseDec = zeros( numIter,numComp );
OptDec = zeros( 0,numComp ); 
CostDesign_CD = []; CostSysFail_CD = [];
numSol_CD = zeros(numIter,1);

for ii = 1:numIter
    baseD_i = RandDec( numDecisions );
    BaseDec(ii,:) = baseD_i;
    [CostDesignDiff_i,CostSysFailProxDiff_i] = evalCostDiff_kN( baseD_i,Cost,CPM_noCost,var,B,varElimOrder );
    OptDec_i = evalOptDecFromProx( CostDesignDiff_i,CostSysFailProxDiff_i,numSubOpt );
    [CostDesign_i,CostSysFail_i] = evalPopKn( OptDec_i,Cost,CPM_noCost,var,B,varElimOrder );
    
    OptDec = [OptDec; OptDec_i]; CostDesign_CD = [CostDesign_CD; CostDesign_i]; CostSysFail_CD = [CostSysFail_CD; CostSysFail_i];
    [OptDec,OptDecUniqueIdx] = unique(OptDec,'rows'); CostDesign_CD = CostDesign_CD(OptDecUniqueIdx); CostSysFail_CD = CostSysFail_CD(OptDecUniqueIdx);
    
    OptDecSortedBoolean = SortNonDominSol_Boolean( CostDesign_CD,CostSysFail_CD );
    OptDec = OptDec(OptDecSortedBoolean,:); CostDesign_CD = CostDesign_CD(OptDecSortedBoolean); CostSysFail_CD = CostSysFail_CD(OptDecSortedBoolean);
   
    plot( CostSysFail_CD,CostDesign_CD,'*' ); grid on; drawnow;
    
end

save(['demoKN_optim_CD_iter' num2str(numIter)])   
%}

%% Draw Results
load demoKN_optim_CD_iter20
load demoKN_optim_gen1500

fontSize_AxisLabel = 17;
fontSize_AxisTick = 14;
fontSize_Legend = 16;

figure;
plot( CostSysFail,CostDesign,'^' )
hold on
plot( CostSysFail_CD,CostDesign_CD,'*' )
grid on

set(gca, 'FontSize', fontSize_AxisTick,'FontName','times new roman')
xlabel( 'System Failure Probability','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )
ylabel( 'Design Cost','Fontsize',fontSize_AxisLabel,'FontName','times new roman' )
% legend( {'Genetic Algorithm' 'Byun and Song (In preparation)'},'Fontsize',fontSize_Legend,'FontName','times new roman','location','northeast')
legend( {'Genetic Algorithm' 'Proxy objective function'},'Fontsize',fontSize_Legend,'FontName','times new roman','location','northeast')

saveas( gcf,'figure/kN_optimResult.emf' )