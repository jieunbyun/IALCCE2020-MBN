function [eventUnspecSet,eventUnspecSetProb] = updateUnspecSet( eventUnspecSet,eventUnspecSetProb,eventIdOut,eventUnspecSetIn,eventUnspecSetInProb )

eventUnspecSet(eventIdOut) = [];
eventUnspecSetProb(eventIdOut) = [];

eventUnspecSet = [eventUnspecSet; eventUnspecSetIn];
eventUnspecSetProb = [eventUnspecSetProb; eventUnspecSetInProb];