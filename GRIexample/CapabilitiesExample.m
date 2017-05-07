%% Load GRI-Mech 3.0 dataset

[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
s = load('GRIMech_modeldata.mat');
modelData = s.GRIMech_modeldata;
dsGRI = startDataset('GRI Mech 3.0');
dsGRI = addData(dsGRI,experimentData,modelData);

%% Check consistency

dsGRI.isConsistent;

%% Plot sensitivity

dsGRI.plotConsistencySensitivity

%% Remove Top 2 

dsGRI.deleteUnit([36,37])

%% Recheck consistency

dsGRI.isConsistent

%% What can explored with a consistent dataset?
%
% * Posterior bounds on Model variables
% * Posterior bounds on dataset QOIs (prediction)
% * Posterior bounds on unmeasured QOIs 
% * Examine distributions and correlations
% * Model optimization


%% Posterior bounds on Model variables
vars = dsGRI.Variables;
nVars = vars.Length;
Opt = B2BDC.Option({'Display',false});
for i1 = 1:2 %nVar
    xi = generateModel(vars.Values(i1));
    dsGRI.setQOI2predict(xi);    
    dsGRI.predictQOI(Opt);
    bndsPosterior(i1) = dsGRI.QOIRange;
end
%% Posterior bounds on dataset QOIs (prediction)

%% Posterior bounds on unmeasured QOIs (prediction)

%% Examine distributions and correlations
xVals = collectSampleDClab(dsGRI);
[h,ax]= plotmatrix(xVals(:,1:5));
[ax.XLim] = deal([-1,1]);
[ax.YLim] = deal([-1,1]);

%% Predict QOIs for xVals
nSamples = size(xVals,1);
yPred = zeros(nSamples,nQOIs);
varList = dsGRI.VarNames;
for i2=1:nQOIs
   qoiModel = dsGRI.DatasetUnits.Values(i2).SurrogateModel;
   yPred(:,i2) = qoiModel.evalActive(xVals,varList);
end
