<<<<<<< Updated upstream
function polyConsisabs(obj,opt)
% POLYCONSISREL(OBJ,OPT) calculate the consistency measure of the dataset
% when all dataset units have polynomial surrogate models and at least one
% of them is not quadratic function.

%  Created: Dec 16, 2016     Wenyu Li

if nargin < 2
   opt = generateOpt;
end
disflag = opt.Display;
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(n_units,1);
bds = obj.calBound;
if opt.AddFitError
   obj.FeasibleFlag = true;
   for i = 1:n_units
      if ~isempty(units(i).SurrogateModel.ErrorStats.absMax)
         abE(i) = units(i).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
bds = bds + [-abE, abE];
db = bds(:,2)-bds(:,1);
allNames = obj.VarNames;
vars = obj.Variables;
ob1 = obj.eval(zeros(1,vars.Length));
tmpc = bds - repmat(ob1',1,2);
cm = generateVar({'CM'},[-max([tmpc(:,1); -tmpc(:,2)]),max(db)],0);
vars = vars.addList(cm);
xbd = vars.calBound;
nall = vars.Length;
tS = zeros(1,nall);
tS(end) = 1;
tPoly = generateModel(-1,tS,vars);
objPoly = tPoly.createSparsePOP;
ineqPolySys = cell(1,2*n_units);
d = 0;
for i = 1:n_units
   if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
      tModel = units(i).SurrogateModel.convertToPoly;
   else
      tModel = units(i).SurrogateModel;
   end
   d = max(d,tModel.Degree);
   tName = tModel.VarNames;
   nM = length(tModel.Coefficient);
   [~,~,id] = intersect(tName,allNames,'stable');
   ts = zeros(nM+1,nall);
   ts(1:end-1,id) = tModel.SupportMatrix;
   ts(end,end) = 1;
   c1 = tModel.Coefficient;
   c1(end+1) = -1;
   p1 = generateModel(c1,ts,vars);
   c2 = -c1;
   c2(end) = -1;
   p2 = generateModel(c2,ts,vars);
   p1 = p1.addConstant(-bds(i,1));
   p2 = p2.addConstant(bds(i,2));
   ineqPolySys{2*i-1} = p1.createSparsePOP;
   ineqPolySys{2*i} = p2.createSparsePOP;
end
param = opt.POPOption.Value;
if param.relaxOrder == -1
   param.relaxOrder = ceil(0.5*d);
end

if disflag
   disp('=======================================================');
   disp('Calculating consistency measure...');
   disp('=======================================================');
end

[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
   sparsePOP(objPoly,ineqPolySys,xbd(:,1),xbd(:,2),param);
x0 = POP.xVectL(1:end-1);
if obj.isFeasiblePoint(x0)
   obj.FeasiblePoint = x0;
   obj.ConsistencyMeasure = [-POP.objValueL -SDPobjValue];
   if disflag
      disp(' ')
      disp('The calculation is done')
      disp(['Consistency LB: ' num2str(-POP.objValueL)])
      disp(['Consistency UB: ' num2str(-SDPobjValue)])
   end
else
   obj.FeasibleFlag = false;
   disp('The calculation is failed, try to use a different POP solver')
end
=======
function polyConsisabs(obj,opt)
% POLYCONSISREL(OBJ,OPT) calculate the consistency measure of the dataset
% when all dataset units have polynomial surrogate models and at least one
% of them is not quadratic function.

%  Created: Dec 16, 2016     Wenyu Li

if nargin < 2
   opt = generateOpt;
end
disflag = opt.Display;
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(n_units,1);
bds = obj.calBound;
if opt.AddFitError
   obj.FeasibleFlag = true;
   for i = 1:n_units
      if ~isempty(units(i).SurrogateModel.ErrorStats.absMax)
         abE(i) = units(i).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
bds = bds + [-abE, abE];
db = bds(:,2)-bds(:,1);
allNames = obj.VarNames;
vars = obj.Variables;
ob1 = obj.eval(zeros(1,vars.Length));
tmpc = bds - repmat(ob1',1,2);
cm = generateVar({'CM'},[-max([tmpc(:,1); -tmpc(:,2)]),max(db)],0);
vars = vars.addList(cm);
xbd = vars.calBound;
nall = vars.Length;
tS = zeros(1,nall);
tS(end) = 1;
tPoly = generateModel(-1,tS,vars);
objPoly = tPoly.createSparsePOP;
ineqPolySys = cell(1,2*n_units);
d = 0;
for i = 1:n_units
   if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
      tModel = units(i).SurrogateModel.convertToPoly;
   else
      tModel = units(i).SurrogateModel;
   end
   d = max(d,tModel.Degree);
   tName = tModel.VarNames;
   nM = length(tModel.Coefficient);
   [~,~,id] = intersect(tName,allNames,'stable');
   ts = zeros(nM+1,nall);
   ts(1:end-1,id) = tModel.SupportMatrix;
   ts(end,end) = 1;
   c1 = tModel.Coefficient;
   c1(end+1) = -1;
   p1 = generateModel(c1,ts,vars);
   c2 = -c1;
   c2(end) = -1;
   p2 = generateModel(c2,ts,vars);
   p1 = p1.addConstant(-bds(i,1));
   p2 = p2.addConstant(bds(i,2));
   ineqPolySys{2*i-1} = p1.createSparsePOP;
   ineqPolySys{2*i} = p2.createSparsePOP;
end
param = opt.POPOption.Value;
if param.relaxOrder == -1
   param.relaxOrder = ceil(0.5*d);
end

if disflag
   disp('=======================================================');
   disp('Calculating consistency measure...');
   disp('=======================================================');
end

[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
   sparsePOP(objPoly,ineqPolySys,xbd(:,1),xbd(:,2),param);
x0 = POP.xVectL(1:end-1);
if obj.isFeasiblePoint(x0)
   obj.FeasiblePoint = x0;
   obj.ConsistencyMeasure = [-POP.objValueL -SDPobjValue];
   if disflag
      disp(' ')
      disp('The calculation is done')
      disp(['Consistency LB: ' num2str(-POP.objValueL)])
      disp(['Consistency UB: ' num2str(-SDPobjValue)])
   end
else
   obj.FeasibleFlag = false;
   disp('The calculation is failed, try to use a different POP solver')
end
>>>>>>> Stashed changes
