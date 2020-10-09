function [CPM, var, B, State] = quantCPM_SF(numLink,linkCoor,EQCoor,linkCapa,numState_M,numState_L,numState_S,state_S,stateDrop_S,probVec_I)

import SF.*
import mbn.*
numCapa = length(linkCapa);

var_ = 0;
% M
var_ = var_+1;
var.M = var_;
State_M = linspace( 6,8.5,numState_M+1 )';
probVec_M = diff( 1-exp( -.76.*State_M ) ); probVec_M = exp( log(probVec_M)-log(sum(probVec_M)) );
CPM{var_} = cpmjoint( var_,(1:numState_M)',probVec_M );
State_M= .5*(State_M(1:end-1)+State_M(2:end));
State{var_}  = State_M;
B{var_} = eye( numState_M );

% L
var_ = var_+1;
var.L = var_;
probVec_L = ones( numState_L, 1) / numState_L;
CPM{var_,1} = cpmjoint( var_,(1:numState_L)',probVec_L );
State{var_,1} = (1:numState_L)';
B{var_,1} = eye( numState_L );

% Sn
state_S = repmat( state_S,numState_S,1 ) - repmat( stateDrop_S*(0:(numState_S-1))',1,length(state_S) );
eventMat_S = (1:numState_S)';
probVec_S = [.7 .2 .1]';
for cc = 1:numLink
    var_ = var_+1;
    var.S(cc,1) = var_;
    CPM{var_} = cpmjoint( var_,eventMat_S,probVec_S );
    State{var_} = state_S;
    B{var_} = eye( numState_S );
end

% In
eventMat_I = [repmat((1:numState_S)',numState_S,1) repelem( (1:numState_S)',numState_S,1 )];
for cc = 1:numLink
    var_ = var_+1;
    var.I(cc,1) = var_;
    CPM{var_} = cpmcond( [var_ var.S(cc)],1,eventMat_I,probVec_I );
    State{var_} = {'Intact' 'Slight Damage' 'Severe Damage'}';
    B{var_} = eye( numState_S );
end

% Xn
eventMat_SML = [repmat((1:numState_M)',numState_L,1) repelem( (1:numState_L)',numState_M,1 )];
eventMat_SML = [repmat((1:numState_S)',numState_M*numState_L,1) repelem( eventMat_SML,numState_S,1 )];
numState_SML = numState_S*numState_M*numState_L;
eventMat_X = [repmat((1:numCapa)',numState_SML,1) repelem( eventMat_SML,numCapa,1 )];
for cc = 1:numLink
    
   var_ = var_+1;
   var.X(cc,1) = var_;
   
   linkCoor_c = linkCoor(cc,:);
   probVec_X_c = [];
   for sml = 1:numState_SML
       stateML_sml = eventMat_SML(sml,:);
       probVec_X_c_sml = PGA2Pf( linkCoor_c,EQCoor(stateML_sml(3),:),State_M(stateML_sml(2)),log( state_S( stateML_sml(1),: ) ) );
       probVec_X_c = [probVec_X_c; probVec_X_c_sml];
   end
   CPM{var_} = cpmcond( [var_ var.S(cc) var.M var.L],1,eventMat_X,probVec_X_c );
   State{var_} = linkCapa(end:-1:1)';
   B{var_} = [eye( numCapa ); 1 1 1; 1 1 0; 0 1 1];  
   
end

% System
var_ = var_+1;
var.System = var_;
CPM{var_} = cpmcond( [var_ var.X'],1 );
State{var_} = {'Survival' 'Failure'}';
B{var_} = eye( 2 );  