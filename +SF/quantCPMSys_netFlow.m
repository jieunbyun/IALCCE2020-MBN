function [eventMat_System,probVec_System] = quantCPMSys_netFlow( eventSurvSet,eventFailSet )

numLink = size( eventSurvSet{1},1 );
numRule_SystemSurv = length( eventSurvSet ); numRule_SystemFail = length( eventFailSet ); 
numRule_System =  numRule_SystemSurv + numRule_SystemFail;
eventMat_System = zeros( numRule_System,1+numLink ); probVec_System = ones( numRule_System,1 );
for ss = 1:numRule_SystemSurv
    rule_s = zeros(1,numLink);
    set_s = eventSurvSet{ss};
    for cc = 1:numLink
        if set_s(cc,1) == set_s(cc,2)
            rule_s(cc) = 4-set_s(cc,1);
        elseif set_s(cc,1)+1 < set_s(cc,2)
            rule_s(cc) = 4;
        else
            if set_s(cc,1) == 2
                rule_s(cc) = 5;
            else
                rule_s(cc) = 6;
            end
        end
    end
    eventMat_System(ss,:) = [1 rule_s];
end

for ss = 1:numRule_SystemFail
    rule_s = zeros(1,numLink);
    set_s = eventFailSet{ss};
    for cc = 1:numLink
        if set_s(cc,1) == set_s(cc,2)
            rule_s(cc) = 4-set_s(cc,1);
        elseif set_s(cc,1)+1 < set_s(cc,2)
            rule_s(cc) = 4;
        else
            if set_s(cc,1) == 2
                rule_s(cc) = 5;
            else
                rule_s(cc) = 6;
            end
        end
    end
    eventMat_System(numRule_SystemSurv+ss,:) = [2 rule_s];
end