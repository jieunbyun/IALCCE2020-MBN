function px = PGA2Pf( aloc,epi,M,lnpga_c )
%{
Compute fail prob. from given EQ scenario
Input:
    aloc: Nn x 2 array of arc-center's location (x,y)
    epi: 1 x 2 array of epicenter's location (x,y)
    M: Scalar of magnitude (6.0-8.5)
    lnpga_c: 1 x (Ns-1) array of log(Sa) 
Output:
    px: 1 x (Ns-1) array of P(X)
%}

loc_v = aloc - epi;
r = sqrt( sum( loc_v.^2,2 ) );
lnpga = -3.512+.904*M-1.328*.5*log(r.^2+(.149*exp(.647*M))^2);

if exp(lnpga) < .068
    sig = .55;
elseif exp(lnpga) < .21
    sig = 0.173-.140*lnpga;
else
    sig = .39;
end

px = normcdf(lnpga_c,lnpga,sig,'upper');
px = diff([0 px]);
px = [px 1-sum(px)];
px = px(end:-1:1);
px = px';
