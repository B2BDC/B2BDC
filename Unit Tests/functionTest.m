% Test script for all functions in B2B package (not including polynomial
% and DClab part)

%  Created: Jan 14, 2016     Wenyu Li

%% generateVar
vName = cell(10,1);
for i = 1:10
   vName{i} = ['v' num2str(i)];
end
H1 = repmat([-1.5,1],10,1);
H2 = repmat([-1,1.5],10,1);
val = zeros(10,1);
v1 = generateVar(vName,H1,val);
v2 = generateVar(vName,H2);

%% functions of VariableList object
% calBound
b1 = v1.calBound;

% addList
v3 = v1.addList(v2);

% deleteVariable
v1 = v1.deleteVariable('v1');
v1 = v1.deleteVariable(3);

% length
length(v1)

% numel
numel(v3)

% makeLHSsample
xSample = v3.makeLHSsample(30);
all(xSample(:) >= -1)
all(xSample(:) <= 1)

% scale
v4 = v3.scale([0.5,2]);
r = rand(10,1);
v5 = v3.scale([0.5*r , 2*r]);

% subDivide
[v6, v7] = v3.subDivide(3, 0.3);
[v8, v9] = v3.subDivide('v3',0.3);

%% generateOpt
opt1 = generateOpt;
opt2 = generateOpt('Display',false, 'ExtraLinFraction',20);

%% generateModel
% a)modelObj = generateModel(GramMatrix, varList)
vName = cell(10,1);
for i = 1:10
   vName{i} = ['v' num2str(i)];
end
H1 = repmat([-1,1],10,1);
v1 = generateVar(vName,H1);
GramMatrix = rand(11);
GramMatrix = GramMatrix + GramMatrix';
m1 = generateModel(GramMatrix, v1);

% eval
xEval = rand(5,10);
y1 = [ones(5,1), xEval] * GramMatrix * [ones(5,1), xEval]';
y1 = diag(y1);
y2 = m1.eval(xEval);
max(abs(y1-y2))

% b)modelObj = generateQModel(vecCoef, varList)
vecCoef = rand(11,1);
m2 = generateModel(vecCoef,v1);

% eval
xEval = rand(5,10);
y1 = [ones(5,1), xEval] * vecCoef;
y2 = m2.eval(xEval);
max(abs(y1-y2))

% c)modelObj = generateQModel(GramN, GramD, varList,k)
GramN = rand(11);
GramN = GramN + GramN';
GramD = zeros(11);
GramD(1,1) = 3;
m3 = generateModel(GramN, GramD, v1, 2);

% eval
xEval = rand(5,10);
y1 = [ones(5,1), xEval] * GramN * [ones(5,1), xEval]';
y1 = diag(y1)/3;
y2 = m3.eval(xEval);
max(abs(y1-y2))

% Dtest  !!!
GramD(1,1) = -2;
m3 = generateModel(GramN, GramD, v1);
GramD(1,1) = 0.3;
m3 = generateModel(GramN, GramD, v1);

% d)modelObj = generateQModel(QN, QD, varList,k)
GramN = rand(11);
GramN = GramN + GramN';
mN = generateModel(GramN, v1);
GramD = zeros(11);
GramD(1,1) = 2;
mD = generateModel(GramD, v1);
m3 = generateModel(mN, mD, v1, 3);
% eval
xEval = rand(5,10);
y1 = [ones(5,1), xEval] * GramN * [ones(5,1), xEval]';
y1 = diag(y1)/2;
y2 = m3.eval(xEval);
max(abs(y1-y2))

%% generateModelbyFit
vName = cell(3,1);
for i = 1:3
   vName{i} = ['v' num2str(i)];
end
H1 = repmat([-2,1],3,1);
v1 = generateVar(vName,H1);
sampleMatrix = rand(4);
sampleMatrix = sampleMatrix + sampleMatrix';
mSample = generateModel(sampleMatrix, v1);
Xsample = rand(100,3);
Ysample = mSample.eval(Xsample)+0.01*rand(100,1);

% modelObj = generateModelbyFit(X, Y, varList)
m1 = generateModelbyFit(Xsample, Ysample, v1);
xTest = rand(50,3);
y1 = mSample.eval(xTest);
y2 = m1.eval(xTest);
max(abs(y1-y2))

% modelObj = generateModelbyFit(X, Y, varList,false,false, q2norm)
m2 = generateModelbyFit(Xsample, Ysample, v1, false, false, 'q2norm');
max(abs(m2.CoefMatrix(:) - sampleMatrix(:)))

% modelObj = generateModelbyFit(X, Y, varList,false,false, rq, K)
m3 = generateModelbyFit(Xsample, Ysample, v1, false, false, 'rq');
m4 = generateModelbyFit(Xsample, Ysample, v1, false, false, 'rq', 12);
y3 = m3.eval(xTest);
y4 = m4.eval(xTest);
max(abs(y1-y3))
max(abs(y1-y4))

%% GenerateDSunit
vName = cell(3,1);
for i = 1:3
   vName{i} = ['v' num2str(i)];
end
H1 = repmat([-2,1],3,1);
v1 = generateVar(vName,H1);
tmpMatrix = rand(4);
tmpMatrix = tmpMatrix + tmpMatrix';
tmpModel = generateModel(tmpMatrix, v1);
unit1 = generateDSunit('test unit', tmpModel, [-1,1]);

% changeBounds
unit2 = unit1.changeBounds([-2,1]);
unit3 = unit1.changeBounds([-2,1,0.2]);

%% generateDataset
dsName = 'Dataset test';
ds1 = generateDataset(dsName);

%% functions of Dataset object
% make up a random but solvable LP problem
vName = cell(8,1);
for i = 1:8
   vName{i} = ['v' num2str(i)];
end
H = repmat([-1,1],8,1);
v = generateVar(vName,H);
xb = v.calBound;
f = rand(8,1);
A = rand(10,8);
b1 = -rand(10,1);
b2 = rand(10,1);
fc = zeros(9,1);
fc(end) = -1;
Ac = [A, ones(10,1);
   -A, ones(10,1)];
[~,cm] = linprog(fc, Ac, [b2;-b1], [],[],[xb(:,1);0], [xb(:,2);1]);
cm = -cm;
[~,pfmin] = linprog(f, [A; -A], [b2;-b1], [],[],xb(:,1), xb(:,2));
[~,pfmax] = linprog(-f, [A; -A], [b2;-b1], [],[],xb(:,1), xb(:,2));
pfmax = -pfmax;

% make up the dataset object
ds = generateDataset('testDS');

% isempty
ds.isempty

% add dataset unit
for i = 1:10
   tmpModel = generateModel([0,A(i,:)], v);
   tmpUnit = generateDSunit(['unit_' num2str(i)],tmpModel,[b1(i), b2(i)]);
   ds.addDSunit(tmpUnit);
end

% calBound
bds = ds.calBound;
bds - [b1, b2]

% calObserve
obs = ds.calObserve;
obs - 0.5*(b1+b2)

% clone
ds2 = ds.clone;

% deleteUnit
dName = {ds2.DatasetUnits.Values.Name};
ds2.deleteUnit(1:3);
ds2.deleteUnit(dName([5,7]));
ds2.DatasetUnits.Values.Name

% changeBounds
ds3 = ds.clone;
n1 = dName([2,8]);
newBD = repmat([-3,3],2,1);
ds3.changeBounds(newBD, n1);
ds3.calBound
ds3.changeBounds(newBD, [4,9]);
ds3.calBound
newBD = rand(ds3.Length,2);
newBD(:,1) = -newBD(:,1);
ds3.changeBounds(newBD);
ds3.calBound - newBD

% findDSunit
unit1 = ds.findDSunit(dName(3));

% isConsistent
opt = generateOpt;
opt.ConsistencyMeasure = 'absolute';
ds.isConsistent(opt);
cm

% isFeasiblePoint
x0 = rand(8,1);
ds.isFeasiblePoint(x0)
x1 = ds.FeasiblePoint;
ds.isFeasiblePoint(x1)

% predictQOI
QOImodel = generateModel([0;f],v);
[QOIrange, xOpt, QOIsens] = ds.predictQOI(QOImodel);
[pfmin pfmax]

% plotQOISensitivity
ds.plotQOISensitivity(QOIsens);

% calVarBoundonF
xBound = ds.calVarBoundonF([]);
ds.Variables.calBound - xBound
xBound2 = ds.calVarBoundonF(2:5);

% calQOIBoundonF
QOIbound = ds.calQOIBoundonF([]);
ds.calBound - QOIbound
QOIbound2 = ds.calQOIBoundonF(dName([3,5,6]));

% plotVarRange
ds.plotVarRange([3,6]);

% greedyDeletionFilter
ds2 = ds.clone;
bds = ds2.calBound;
% make unit 7 or 9 the cause for inconsistency
ds2.changeBounds(0.3*bds(7,:)-2,7);
ds2.changeBounds(0.6*bds(9,:)+4,9);
ds2.isConsistent(opt);
ds3 = ds2.clone;
[tepUnit, dgap] = ds3.greedyDeletionFilter(opt);
tepUnit
dgap

% makeConsistency
ds3 = ds2.clone;
ds3.makeConsistency(20,opt);

% plotMinBoundChange
ds3 = ds2.clone;
ds3.plotMinBoundChange(opt)

% selfInconsisAnalysis
ds3 = ds2.clone;
ds3.changeBounds([-25,-22],3);
dsNew = ds3.selfInconsisAnalysis(opt);