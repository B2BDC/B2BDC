[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
s = load('GRIMech_modeldata.mat');
modelData = s.GRIMech_modeldata;
dsGRI = startDataset('GRI Mech 3.0');
dsGRI = addData(dsGRI,experimentData,modelData);
dsGRI.deleteUnit([37])
dsGRI.isConsistent;

%% Runs = 500; BlockSize = 10000  -> 15 minutes
% Runs = 24*500; BlockSize = 10000  -> 6 hours
xCenters = -0.99:0.02:0.99;
xStats = 2;
YM = [];
yCenters = [];
yStats = [];
N.Runs = 12000;
N.BlockSize = 10000;
N.ReturnSamples = 10000;
[xEmp, yEmp, Xvals] = collectSamplesNew(dsGRI,[],xCenters,xStats,YM,yCenters,yStats,N);
