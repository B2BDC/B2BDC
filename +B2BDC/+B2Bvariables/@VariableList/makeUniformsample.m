function xCand = makeUniformsample(obj,n)

% XCAND = MAKEUNIFORMSAMPLE(OBJ,N) returns a n-by-nVar sample matrix. The method
% generates the sample with rand function.

H = obj.calBound;
nVar = size(H,1);
dx = diff(H,[],2);
xCand = repmat(H(:,1)',n,1)+repmat(dx',n,1).*rand(n,nVar);