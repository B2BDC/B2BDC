% Test consistency measure, prediction with NNmodel (neural network)
% August 21, 2018     Wenyu Li

%% Create a dataset with quadratic or network models
clear;clc;
seedID = 17;
rng(seedID);
nAll = 10;
nAct = 6;
nQ = 10;
vAll = generateVar([],repmat([0 1],nAll,1));
ds1 = generateDataset('Quadratic DS');
ds2 = generateDataset('Network DS');
ndata = 1e3;
for i = 1:nQ
   idV = randperm(nAll,nAct);
   Q = rand(nAct+1);
   Q = 0.5*(Q+Q');
   tmpV = vAll.makeSubset(idV);
   QModel = generateModel(Q,tmpV);
   xdata = tmpV.makeLHSsample(ndata);
   ydata = QModel.eval(xdata);
   NNModel = generateModelbyFit(xdata,ydata,tmpV,'nn',1);
   LB = Q(1,1)-rand(1);
   UB = Q(1,1)+rand(1);
   tmpUnit1 = generateDSunit(['unit' num2str(i)],QModel,[LB UB]);
   tmpUnit2 = generateDSunit(['unit' num2str(i)],NNModel,[LB UB]);
   ds1.addDSunit(tmpUnit1);
   ds2.addDSunit(tmpUnit2);
end
opt = generateOpt('Display',false,'AddFitError',false);

%% Consistency measure
ds1.clearConsis;
ds2.clearConsis;
ds1.isConsistent(opt);
ds2.isConsistent(opt);

%% generating prediction QOI
idV = randperm(nAll,nAct);
Q = rand(nAct+1);
Q = 0.5*(Q+Q');
tmpV = vAll.makeSubset(idV);
Qpred = generateModel(Q,tmpV);
xdata = tmpV.makeLHSsample(ndata);
ydata = Qpred.eval(xdata);
NNpred = generateModelbyFit(xdata,ydata,tmpV,'nn',1);

%% making prediction
[p1,s1] = ds1.predictQOI(Qpred,opt);
[p2,s2] = ds2.predictQOI(NNpred,opt);

%% finite difference gradient
clear;clc
nVar = 6;
xdata = rand(1000,nVar);
ydata = rand(1000,1);
vars = generateVar([],repmat([0 1],nVar,1));
NNModel = generateModelbyFit(xdata,ydata,vars,'nn',1,[8 3 2]);
QModel = generateModelbyFit(xdata,ydata,vars,'q2norm',1);
CoefQ = QModel.CoefMatrix;
fg1 = NNModel.NetGradientFcn;
fg2 = @(x) 2*(CoefQ(2:end,2:end)*x+CoefQ(2:end,1));
x1 = rand(nVar,1);
dydx1 = zeros(nVar,1);
dydx2 = zeros(nVar,1);
dx = 1e-10;
for i = 1:nVar
   x1p = x1; x1p(i) = x1p(i)+dx;
   x1n = x1; x1n(i) = x1n(i)-dx;
   dydx1(i) = (NNModel.eval(x1p')-NNModel.eval(x1n'))/2/dx;
   dydx2(i) = (QModel.eval(x1p')-QModel.eval(x1n'))/2/dx;
end
ddy1 = abs(fg1(x1)-dydx1);
ddy2 = abs(fg2(x1)-dydx2);