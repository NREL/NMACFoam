singlePhaseTransportModel laminarTransport(U, phi);

const int polySize=5;

Polynomial<polySize> muCoeffs = 
laminarTransport.lookup
(
     "muCoeffs<" + Foam::name(polySize) + '>'
);

Polynomial<polySize> rhoCoeffs = 
laminarTransport.lookup
(
     "rhoCoeffs<" + Foam::name(polySize) + '>'
);

Polynomial<polySize> DCoeffs = 
laminarTransport.lookup
(
     "DCoeffs<" + Foam::name(polySize) + '>'
);
