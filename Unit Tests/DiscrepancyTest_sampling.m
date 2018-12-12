clear;clc
opt = generateOpt('Display',false,'AddFitError',false,'LocalStart',5);
opt.ConsistencyMeasure = 'absolute';
% opt.Prediction = 'inner';
%% Feasible point check (multiple correction)
nVar = 6;
actVar = 4;
nQOI = 15;
ds1 = generateRandomDataset(nQOI,nVar,actVar);
sv = rand(nQOI,2);
sname = {'Scenario 1','Scenario 2'};
ds1 = ds1.addScenarioParameter(sv,sname);
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
ds1.setModelDiscrepancy(QOIindex,f1,10);
ds1.setParameterDiscrepancy(varIndex,f2(varIndex),Hv,10);
nX = ds1.Variables.Length+ds1.ModelDiscrepancy.Variables.Length+ds1.ParameterDiscrepancy.Variables.Length;
% points with nVar dimension
nS = 10;
x1 = ds1.Variables.makeLHSsample(nS);
[flag1,newx1] = ds1.isFeasiblePoint(x1);
% points with full dimension
x2 = randn(nS,nX);
x2(:,1:nVar) = ds1.Variables.makeLHSsample(nS);
[flag2,newx2] = ds1.isFeasiblePoint(x2);
% consistency
ds1.isConsistent(opt);
xx = ds1.FeasiblePoint';
ds1.clearConsis;
[flag,xnew] = ds1.isFeasiblePoint(xx);

%% HR samples
ds1.clearConsis;
opt.SampleOption.StepInterval = 50;
nS = 100;
% sample with feasible point
xs1 = ds1.collectSamples(nS,[],opt);
flag1 = ds1.isFeasiblePoint(xs1.x);
sum(flag1)
ds1.clearConsis;
flag2 = ds1.isFeasiblePoint(xs1.x(randperm(nS,1),1:xs1.dimension(1)));
% sample with a possible infeasible point
xStart = randn(1,sum(xs1.dimension));
ds1.collectSamples(1,xStart,opt);

%% Modified HR samples
ds1.clearConsis;
opt.SampleOption.StepInterval = 50;
nS = 100;
% sample with feasible point
xs1 = ds1.collectHRsamples(nS,[],opt);
flag1 = ds1.isFeasiblePoint(xs1.x);
sum(flag1)
ds1.clearConsis;
flag2 = ds1.isFeasiblePoint(xs1.x(randperm(nS,1),1:xs1.dimension(1)));
% sample with a possible infeasible point
xStart = randn(1,sum(xs1.dimension));
ds1.collectHRsamples(1,xStart,opt);