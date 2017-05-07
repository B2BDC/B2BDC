function  [xEmp,yEmp,xVals] = collectSampleDClab(obj,opt)
% Gateway function to collect a large number of samples from the feasible
% set defined by the dataset obj through DClab

% Created: July 20, 2015     Wenyu Li

% filename = which('cvxquadratic2normfit.m');
% pth1 = fileparts(filename);
% pth2 = fullfile(pth1,'DClabV2');
% pth3 = fullfile(pth1,'GRIexample');
% addpath(genpath(pth2));
% addpath(pth3);

units = obj.DatasetUnits.Values;
n = length(units);
RO = cell(n,1);
RM = cell(n,1);
Pair = cell(n,1);
vars = obj.Variables.Values;
xmin = min([vars.LowerBound]);
xmax = max([vars.UpperBound]);
% allNames = {};
% p0 = {};
% for i = 1:n
%    varData = units(i).SurrogateModel;
%    j = length(varData.VarNames);
%    allNames(end+1:end+j,1) = varData.VarNames;
%    varValue = varData.Variables.Values;
%    p0(end+1:end+j,1) = {varValue.NominalValue};
% end
% [allNames,idx] = unique(allNames);
% p0 = p0(idx);

for i = 1:n
   d = units(i).ObservedValue;
   u.value = [units(i).LowerBound-d, units(i).UpperBound-d];
   u.type = 'absolute';
   u.transformation='none';
   
   RO{i} = DClab.ResponseObservation(d,u);
   
   pNames = units(i).SurrogateModel.VarNames;
   pRanges = {};
   for j = 1:length(pNames)
      pRanges(j,1) = {[-inf inf]};
   end
   
   paramDomain = struct('name',pNames,'range',pRanges);
   
   RM{i} = DClab.ResponseModel(units(i).SurrogateModel.CoefMatrix,paramDomain);
   
   RM{i} = set(RM{i},'name',units(i).Name);
   
   name = units(i).Name;
   Pair{i,1} = DClab.ModelAndObservationPair(RO{i},RM{i},name);
end

n = length(vars);
P = cell(n,1);
for i = 1:n
   name = vars(i).Name;
   u.type = 'absolute'; u.transformation = 'none';
   u.value = [vars(i).LowerBound, vars(i).UpperBound]-vars(i).NominalValue;
   P{i} = DClab.FreeParameter(name,vars(i).NominalValue,u);
end

FreeParam = vertcat(P{:});

D = DClab.DCDataset(vertcat(Pair{:}),FreeParam);
D.name = obj.Name;

ctest = DClab.ConsistencyTest(D);

if nargin < 2
   xInit = ctest.LBx;
   xCenters = xmin+0.01:0.01:xmax-0.01;
   xStats = 2;
   N.BlockSize = 10000;
   N.Runs = 30;
   N.ReturnSamples = 3000;
   YM = [];
   yCenters = [];
   yStats = [];
end
[xEmp,yEmp, xVals] = ...
   collectSamples(D,xInit,xCenters,xStats,YM,yCenters,yStats,N,ctest.LB);