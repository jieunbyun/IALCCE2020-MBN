function db = RandDec( Nd )
%{
get decision rule by random generation

INPUT:
    Nd: M x 1 array of decision alternatives #
OUTPUT:
    db: 1 x M array with basis decision
%}

M = length(Nd);
db = zeros(1,M);
for ii = 1:M
    db(ii) = randsample( Nd(ii),1 );
end