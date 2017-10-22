% Loads Friday18December.mat and creates a dataset in the B2BDC_v0.6 format
% Compare to corresponding DCLabV2 dataset.
% Arun Hegde
% To do list
% Convert to function that produces random feasible datasets.
%
load Friday18December;

dsB2BTest = startDataset('Friday18DecemberDataset');
nY = numel(MatsVars);
 
for kk = 1:nY
    modelVar = MatsVars(kk).Vars;
    n = size(modelVar,1);
    varLB = -ones(n,1);
    varUB = ones(n,1);
    obs = zeros(n,1);
    varList = B2BDC.B2Bvariables.VariableList();
    for j = 1:n
        modelVariable = B2BDC.B2Bvariables.ModelVariable(modelVar{j},...
            varLB(j),varUB(j),obs(j));
        varList = varList.add(modelVariable);
    end
    coefMat = MatsVars(kk).Matrix;
    unitModel = B2BDC.B2Bmodels.QModel(coefMat,varList);
    unitName = ['Y' int2str(kk)];
    unitLB = MatsVars(kk).YBounds(1);
    unitUB = MatsVars(kk).YBounds(2);
    unitVal = (unitLB+unitUB)/2;
    dsUnit = B2BDC.B2Bdataset.DatasetUnit(unitName,unitModel,unitLB,unitUB,unitVal);
    dsB2BTest.addDSunit(dsUnit);
end

dsB2BTest.isConsistent
%% Test the predictions of the variables on the feasible set

vars = dsB2BTest.Variables;
nVars = vars.Length;
Opt = B2BDC.Option({'Display',false});
for i1 = 1:nVars
    xi = generateModel(vars.Values(i1));
    dsB2BTest.setQOI2predict(xi);    
    dsB2BTest.predictQOI(Opt);
    bndsPosterior(i1) = dsB2BTest.QOIRange;
end

%% Sampling from the feasible set
xCenters = [-0.99:0.01:0.99];
xStats = 1;
N.BlockSize = 10000;
N.Runs = 30;
N.ReturnSamples = 1000;
[xEmp,yEmp,Xvals] = collectSamplesNew(dsB2BTest,[],xCenters,xStats,[],[],[],N);