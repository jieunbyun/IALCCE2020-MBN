function D = GenMut( Nmu,Nd,D,rmu )
%{
Generate Nmu decision rules by mutation from the ones in D
INPUT: 
    Nmu: scalar of mutated population # / generation
    Nd: 1 x Nx array of decision alternative #
    D: Nsol x Nx array of decision rules
    rmu: scalar of genetic change ratio
OUTPUT:
    D: Nmu x Nx array of decision rules
%}

Nsol = size(D,1); Nx = length(Nd); D2 = zeros(Nmu,Nx);
for ii = 1:Nmu
    tmp_ = randsample(Nsol,1);
    tmp2_ = randsample(2,Nx,'true',[1-rmu rmu])-1;
    d = D(tmp_,:);
    tmp2_ = find(tmp2_);
    for jj = 1:length(tmp2_)
        tj = tmp2_(jj);
        d(tj) = randsample( Nd(tj),1 );
    end
    D2(ii,:) = d;
end
D = D2;