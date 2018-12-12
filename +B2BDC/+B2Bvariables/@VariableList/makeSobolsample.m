function xCand = makeSobolsample(obj,n)

% XCAND = SOBOLSAMPLE(OBJ,N) returns a n-by-nVar sample matrix. The method
% generates the sample with sobol sequence.

nVar = obj.Length;
sob = sobolset(nVar,'Skip',1e3*randi(1e3),'Leap',1e2*randi(5e2));
sob = scramble(sob,'MatousekAffineOwen');
xx = net(sob,n);
bds = obj.calBound;
dx = bds(:,2) - bds(:,1);
xCand = repmat(bds(:,1)',n,1) + xx.*repmat(dx',n,1);