function [QOIrange, QOISensitivity, xOpt] = preQOIpoly(obj,QOIobj,opt,rflag)
% RESULT = PREQOIPOLY(OBJ,OPT) calculates the QOI prediction for a
% dataset with all polynomial surrogate models and at least one of them is
% not quadratic.

%  Created: Dec 23, 2016     Wenyu Li

if nargin < 3
   opt = generateOpt;
end
if nargin < 4
   rflag = false;
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(n_units,1);
bds = obj.calBound;
if opt.AddFitError
   for i = 1:n_units
      if ~isempty(units(i).SurrogateModel.ErrorStats.absMax)
         abE(i) = units(i).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
bds = bds + [-abE, abE];
allNames = obj.VarNames;
vars = obj.Variables;
nall = vars.Length;
xbd = vars.calBound;
if isa(QOIobj,'B2BDC.B2Bmodels.QModel')
   objModel = QOIobj.convertToPoly;
else
   objModel = QOIobj;
end
objPoly = objModel.createSparsePOP;
ineqPolySys = cell(1,2*n_units);
d = 0;
if rflag
   ncount = n_units-1;
else
   ncount = n_units;
end
for i = 1:ncount
   if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
      tModel = units(i).SurrogateModel.convertToPoly;
   else
      tModel = units(i).SurrogateModel;
   end
   d = max(d,tModel.Degree);
   tName = tModel.VarNames;
   nM = length(tModel.Coefficient);
   [~,~,id] = intersect(tName,allNames,'stable');
   ts = zeros(nM,nall);
   ts(:,id) = tModel.SupportMatrix;
   c1 = tModel.Coefficient;
   p1 = generateModel(c1,ts,vars);
   c2 = -c1;
   p2 = generateModel(c2,ts,vars);
   p1 = p1.addConstant(-bds(i,1));
   p2 = p2.addConstant(bds(i,2));
   ineqPolySys{2*i-1} = p1.createSparsePOP;
   ineqPolySys{2*i} = p2.createSparsePOP;
end
param = opt.POPOption.Value;
if param.relaxOrder == -1
   param.relaxOrder = ceil(0.5*d);O
end

[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
   sparsePOP(objPoly,ineqPolySys,xbd(:,1),xbd(:,2),param);

QOIrange.min = [SDPobjValue POP.objValueL];
xOpt = POP.xVectL;

objPoly.coef = -objPoly.coef;
[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
   sparsePOP(objPoly,ineqPolySys,xbd(:,1),xbd(:,2),param);
QOIrange.max = [-POP.objValueL -SDPobjValue];
xOpt = [xOpt, POP.xVectL];

QOISensitivity = [];