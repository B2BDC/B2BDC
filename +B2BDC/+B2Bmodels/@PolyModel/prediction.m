function prediction(obj,opt)
% To predict the a B2BDC.B2Bmodels.Polymodel object through sparsePOP
% package

% Created: Nov 27, 2015     Wenyu Li

if nargin < 2
   opt = generateOpt;
end
objPoly.typeCone = 1;
pBasis = obj.Basis;
objPoly.degree = pBasis.calDegree;
objPoly.dimVar = pBasis.Dimension;
objPoly.noTerms = pBasis.Length;
objPoly.supports = pBasis.Value;
objPoly.coef = obj.Coefficient;
ineqPolySys = {};
varbound = obj.Variables.calBound;
lbd = (varbound(:,1))';
ubd = (varbound(:,2))';
param.relaxOrder = opt.SOSrelaxOrder;
param.printFileName = 0;
% param.sparseSW = 0;
[param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] =...
   sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
clc
obj.PredictionRange.min = [SDPobjValue, POP.objValue];
objPoly.coef = -obj.Coefficient;
[param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] =...
   sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
clc
obj.PredictionRange.max = [-POP.objValue, -SDPobjValue];
if opt.Display
   disp(['The minimum is within [' num2str(obj.PredictionRange.min(1)) '  ' num2str(obj.PredictionRange.min(2)) ']'])
   disp(['The maximum is within [' num2str(obj.PredictionRange.max(1)) '  ' num2str(obj.PredictionRange.max(2)) ']'])
end
