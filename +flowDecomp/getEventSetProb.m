function eventSetProb = getEventSetProb( eventSet,linkProb )

if ~isempty(eventSet)
    numLink = size(eventSet{1},1);
end
   
numEventSet = length(eventSet);
eventSetProb = zeros(numEventSet,1);
for ee = 1:numEventSet
    event = eventSet{ee};
    prob_e = log(1);
    for ll = 1:numLink
        linkBnd_l = event(ll,:);
        prob_e = prob_e + log( sum(linkProb(ll,linkBnd_l(1):linkBnd_l(2)) ) );
    end
    prob_e = exp(prob_e);
    eventSetProb(ee) = prob_e;
end