function polyConsisrel(obj,opt)

% This function is used to calculate the relative consistency measure for a
% dataset consists of only polynomial surrogate and prediction models.

%  Created: Nov 30, 2015     Wenyu Li

vars = obj.Variables;
optimalVar = generateVar({'eps'},[-1e10,1e10]);
newVar = vars.addList(optimalVar);
nVar = vars.Length;
nUnit = obj.Length;
ineqPoly = cell(2*nUnit,1);
bds = obj.calBound;
ob = obj.calObserve;
unc = bds-repmat(ob,1,2);
tep_eps = zeros(1,nVar+1);
tep_eps(end) = 1;
for i = 1:nUnit
   pm1 = obj.DatasetUnits.Values(i).SurrogateModel;
   pm2 = pm1.clone;
   pm1 = pm1.expandDimension(newVar);
   pm2 = pm2.expandDimension(newVar);
   pm2.Coefficient = -pm2.Coefficient;
   [~,id] = intersect(pm1.Basis.Value,zeros(1,nVar+1),'rows');
   pm1.Coefficient(id) = pm1.Coefficient(id)-bds(i,1);
   [~,id] = intersect(pm2.Basis.Value,zeros(1,nVar+1),'rows');
   pm2.Coefficient(id) = pm2.Coefficient(id)+bds(i,2);
   pm1.Basis.Value = [pm1.Basis.Value; tep_eps];
   pm1.Coefficient(end+1) = unc(i,1);
   pm2.Basis.Value = [pm2.Basis.Value; tep_eps];
   pm2.Coefficient(end+1) = -unc(i,2);
   ineqPoly{2*i-1}.typeCone = 1;
   ineqPoly{2*i-1}.degree = pm1.Basis.calDegree;
   ineqPoly{2*i-1}.dimVar = nVar+1;
   ineqPoly{2*i-1}.noTerms = pm1.Basis.Length;
   ineqPoly{2*i-1}.supports = pm1.Basis.Value;
   ineqPoly{2*i-1}.coef = pm1.Coefficient;
   ineqPoly{2*i}.typeCone = 1;
   ineqPoly{2*i}.degree = pm2.Basis.calDegree;
   ineqPoly{2*i}.dimVar = nVar+1;
   ineqPoly{2*i}.noTerms = pm2.Basis.Length;
   ineqPoly{2*i}.supports = pm2.Basis.Value;
   ineqPoly{2*i}.coef = pm2.Coefficient;
end
objPoly.typeCone = 1;
objPoly.degree = 1;
objPoly.dimVar = nVar+1;
objPoly.noTerms = 1;
objPoly.supports = tep_eps;
objPoly.coef = -1;
xbds = newVar.calBound;
param.relaxOrder = opt.SOSrelaxOrder;
[param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] =...
   sparsePOP(objPoly,ineqPoly,xbds(:,1),xbds(:,2),param);
clc;

