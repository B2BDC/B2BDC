clear;clc
opt = generateOpt('Display',false,'AddFitError',false,'LocalStart',5);
opt.ConsistencyMeasure = 'absolute';
% opt.Prediction = 'inner';
%% Parameter discrepancy with model discrepancy (multiple correction)
nVar = 6;
actVar = 4;
nQOI = 16;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
QQ.model = ds1.DatasetUnits.Values(nQOI).SurrogateModel;
ds1.deleteUnit(nQOI);
sv = rand(nQOI,2);
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv(1:end-1,:),sname);
f1{1} = @(x) [ones(size(x,1),1) x(:,1) exp(x(:,2))];
f1{2} = @(x) [0.5*x(:,1).*x(:,2) -exp(x(:,1))];
QOIindex{1} = 1:8;
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
xtest = randn(1,16);
save xtest xtest
QQ.group = 1;
QQ.sv = sv(end,:);
Qcorrect.GroupIndex = QQ.group;
Qcorrect.Value = sv(end,:);
[x,y] = calculatePrediction(ds1,QQ,QOIindex,f1,f2,sv(1:end-1,:),10,Hv);
ds1.setModelDiscrepancy(QOIindex,f1,10);
ds1.setParameterDiscrepancy(varIndex,f2(varIndex),Hv,10);
[yy,ss,xx] = ds1.predictQOI(QQ.model,opt,Qcorrect);
