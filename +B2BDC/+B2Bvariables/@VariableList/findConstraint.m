function [varList,constraint] = findConstraint(obj,xg,xb,np)
% VARLIST = FINDLEVEL2CONSTRAINT(OBJ,XG,XB,NC) returns a new variableList
% object with added linear constraints. The number of constraints is
% specified by nc.

%  Created: Oct 11, 2017     Wenyu Li

cuttingDegree = 0.1;
[nS,nVar] = size(xg);
xc = xg - repmat(mean(xg),nS,1);
[V,D] = eig(xc' * xc);
d = diag(D);
d = d/sum(d);
[~,id] = sort(d,'descend');
V = V(:,id);
A = zeros(np,nVar);
LB = zeros(np,1);
UB = zeros(np,1);
for i = 1:np
   tv = V(:,i);
   tv = tv/max(abs(tv));
   A(i,:) = tv;
   td = xb*tv;
   td_p = td(td>=0);
   td_n = td(td<0);
   LB(i) = quantile(td_n,1-cuttingDegree);
   UB(i) = quantile(td_p,cuttingDegree);
end
varList= obj.addLinearConstraint(A,LB,UB);
constraint.A = A;
constraint.LB = LB;
constraint.UB = UB;


