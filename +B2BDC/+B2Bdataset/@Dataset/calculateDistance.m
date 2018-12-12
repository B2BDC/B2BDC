function [d,xmin] = calculateDistance(obj,x,opt,Curve)
% D = CALCULATEDISTANCE(OBJ,X,OPT,CURVE) calculates the distance of the given point x
% from the feasible set, defined as minimum sqrt((x'-x)'(x'-x)) over x' in
% F. It returns infinity if the dataset is inconsistent.

%  Created: Sep 16, 2018     Wenyu Li

if size(x,1) == 1
   x = x';
end
nVar = obj.Variables.Length;
if nVar ~= length(x)
   error('The given point has a wrong dimension!');
end
if nargin < 3
   opt = generateOpt('Display',false,'Prediction','inner');
end
if ~obj.isConsistent(opt)
   d = inf;
elseif nargin > 3
   [d,xmin] = calculateDistanceCurve(obj,x,opt,Curve);
elseif ~obj.ParameterDiscrepancyFlag
   Coef = zeros(nVar+1);
   Coef(2:end,2:end) = eye(nVar);
   Coef(1,2:end) = x;
   Coef(2:end,1) = x;
   Coef(1,1) = x'*x;
   Qpred = generateModel(Coef,obj.Variables);
   [dd,~,xopt] = obj.predictQOI(Qpred,opt);
   xmin = xopt.min;
   d = sqrt(dd.min); 
else
   sv = obj.getScenarioParameter;
   svmin = min(sv);
   dsv = max(sv)-min(sv);
   C.Value = repmat(svmin,101,1)+repmat(dsv,101,1).*repmat((0:0.01:1)',1,length(svmin));
   C.function = @(x) zeros(101,length(x));
   [d,xmin] = calculateDistanceCurve(obj,x,opt,C);
%    error('Scenario value needs to be specified for parameter discrepancy');
end