function xCand = makeSobolsample(obj,n)

% XCAND = SOBOLSAMPLE(OBJ,N) returns a n-by-nVar sample matrix. The method
% generates the sample with sobol sequence.

nVar = obj.Length;
sob = sobolset(nVar,'Skip',1e3*ceil(100*rand(1)),'Leap',1e2*ceil(100*rand(1)));
sob = scramble(sob,'MatousekAffineOwen');
xx = net(sob,n);
bds = obj.calBound;
dx = bds(:,2) - bds(:,1);
xCand = repmat(bds(:,1)',n,1) + xx.*repmat(dx',n,1);