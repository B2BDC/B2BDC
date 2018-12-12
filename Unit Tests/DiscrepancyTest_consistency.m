clear;clc
opt = generateOpt('Display',false,'AddFitError',false,'LocalStart',5);
opt.ConsistencyMeasure = 'absolute';
%% Model discrepancy test (single correction)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(15,2);
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
f1 = @(x) [ones(size(x,1),1) x(:,1) exp(x(:,2))];
% validation calculation
xtest = rand(1,10);
save xtest xtest
[x,y,lambda] = calculateConsistency(ds1,[],f1,[],sv,10);
ds1.setModelDiscrepancy([],f1,10);
ds1.isConsistent(opt);
xcom = [x(ds1.Variables.Length+1:end) ds1.ModelDiscrepancy.FeasiblePoint];
ycom = [y ds1.ConsistencyMeasure(1)];
load y1
load y2
dy = y1-y2;
%% Model discrepancy test (multiple corrections)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(15,2); 
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
QOIgroup{1} = 1:5;
QOIgroup{2} = 7:13;
f1{1} = @(x) [ones(size(x,1),1) x(:,1) exp(x(:,2))];
f1{2} = @(x) [ones(size(x,1),1) x(:,1).^0.8 x(:,1)-x(:,2)];
% validation calculation
xtest = rand(1,13);
save xtest xtest
[x,y,lambda] = calculateConsistency(ds1,QOIgroup,f1,[],sv,10);
ds1.setModelDiscrepancy(QOIgroup,f1,10);
ds1.isConsistent(opt);
xcom = [x(ds1.Variables.Length+1:end) ds1.ModelDiscrepancy.FeasiblePoint];
ycom = [y ds1.ConsistencyMeasure(1)];
load y1
load y2
dy = y1-y2;
%% Parameter discrepancy test (no QOI correction)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(15,2); 
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
f2 = cell(nVar,1);
f2{2} = @(x) [x(:,1).^2, x(:,2)];
f2{5} = @(x) [x(:,1).*x(:,2)];
varIndex = [5 2];
oldH = ds1.Variables.calBound;
Hv = repmat([-inf inf],ds1.Variables.Length,1);
Hv(varIndex,:) = 2*oldH(varIndex,:);
% validation calculation
xtest = rand(1,10);
save xtest xtest
[x,y,lambda] = calculateConsistency(ds1,[],[],f2,sv,10,Hv);
ds1.setParameterDiscrepancy(varIndex,f2(varIndex),Hv,10);
ds1.isConsistent(opt);
xcom = [x(ds1.Variables.Length+1:end) ds1.ParameterDiscrepancy.FeasiblePoint];
ycom = [y ds1.ConsistencyMeasure(1)];
load y1
load y2
dy = y1-y2;
%% Parameter discrepancy with model discrepancy (single correction)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(15,2); 
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
f1 = @(x) [ones(size(x,1),1) x(:,1) exp(x(:,2))];
f2 = cell(nVar,1);
f2{2} = @(x) [x(:,1).^2, x(:,2)];
f2{3} = @(x) [0.5*x(:,2) 0.2*x(:,1)];
f2{5} = @(x) [x(:,1).*x(:,2)];
varIndex = [5 2 3];
oldH = ds1.Variables.calBound;
Hv = repmat([-inf inf],ds1.Variables.Length,1);
Hv(varIndex,:) = 2*oldH(varIndex,:);
% validation calculation
xtest = rand(1,15);
save xtest xtest
[x,y,lambda] = calculateConsistency(ds1,[],f1,f2,sv,10,Hv);
ds1.setModelDiscrepancy([],f1,10);
ds1.setParameterDiscrepancy(varIndex,f2(varIndex),Hv,10);
ds1.isConsistent(opt);
xcom = [x(ds1.Variables.Length+1:end) [ds1.ModelDiscrepancy.FeasiblePoint;  ds1.ParameterDiscrepancy.FeasiblePoint]];
ycom = [y ds1.ConsistencyMeasure(1)];
load y1
load y2
dy = y1-y2;
%% Parameter discrepancy with model discrepancy (multiple correction)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(15,2); 
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
f1{1} = @(x) [ones(size(x,1),1) x(:,1) exp(x(:,2))];
f1{2} = @(x) [0.5*x(:,1).*x(:,2) -exp(x(:,1))];
QOIindex{1} = 1:10;
QOIindex{2} = 13:15;
f2 = cell(nVar,1);
f2{2} = @(x) [x(:,1).^2, x(:,2)];
f2{3} = @(x) [0.5*x(:,2) 0.2*x(:,1)];
f2{5} = @(x) [x(:,1).*x(:,2)];
varIndex = [5 2 3];
oldH = ds1.Variables.calBound;
Hv = repmat([-inf inf],ds1.Variables.Length,1);
Hv(varIndex,:) = 2*oldH(varIndex,:);
% validation calculation
xtest = rand(1,17);
save xtest xtest
[x,y,lambda] = calculateConsistency(ds1,QOIindex,f1,f2,sv,10,Hv);
ds1.setModelDiscrepancy(QOIindex,f1,10);
ds1.setParameterDiscrepancy(varIndex,f2(varIndex),Hv,10);
ds1.isConsistent(opt);
xcom = [x(ds1.Variables.Length+1:end) [ds1.ModelDiscrepancy.FeasiblePoint;  ds1.ParameterDiscrepancy.FeasiblePoint]];
ycom = [y ds1.ConsistencyMeasure(1)];
% load y1
% load y2
% dy = y1-y2;
%% Consistency Outer Bound
clear;clc
nVar = 12;
actVar = 8;
nQOI = 30;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(nQOI,2); 
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
opt = generateOpt('Display',false,'AddFitError',false,'LocalStart',5);
f1 = @(x) ones(size(x,1),1);
ds1.setModelDiscrepancy([],f1,10);
ds1.isConsistent(opt);
basis = f1(sv);
vExtra = ds1.ModelDiscrepancy.Variables;
nExtra = vExtra.Length;
ds2 = generateDataset('Test');
bds = ds1.calBound;
obs = ds1.calObserve;
for i = 1:nQOI
   v0 = ds1.DatasetUnits.Values(i).SurrogateModel.Variables;
   v0 = v0.addList(vExtra);
   Coef = ds1.DatasetUnits.Values(i).SurrogateModel.CoefMatrix;
   newCoef = zeros(size(Coef,1)+nExtra);
   newCoef(1:end-nExtra,1:end-nExtra) = Coef;
   newCoef(end-nExtra+1:end,1) = 0.5*basis(i,:);
   newCoef(1,end-nExtra+1:end) = 0.5*basis(i,:);
   newModel = generateModel(newCoef,v0);
   tUnit = generateDSunit(ds1.DatasetUnits.Values(i).Name,newModel,bds(i,:),obs(i));
   ds2.addDSunit(tUnit);
end
ds2.isConsistent(opt);