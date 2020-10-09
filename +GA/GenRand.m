function D = GenRand( Npop,Nd )
%{
Generate Npop decision rules randomly
INPUT: 
    Npop: scalar of population # / generation
    Nd: 1 x Nx array of decision alternative #
OUTPUT:
    D: Npop x Nx matrix with generated rules
%}

Nx = length(Nd);
D = zeros(Npop,Nx);
for ii = 1:Npop
    d_ = zeros( 1,Nx );
    for jj = 1:Nx
        d_(jj) = randsample(Nd(jj),1);
    end
    D(ii,:) = d_;
end
