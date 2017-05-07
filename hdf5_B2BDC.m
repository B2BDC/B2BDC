function dsNew = hdf5_B2BDC(lc)
% Gateway function used to transform a hdf5 file into a 
% B2BDC.B2Bdataset.Dataset object.
% lc - The location of the input hdf5 file

% Created: Oct 13, 2015     Wenyu Li

%% Load dataset hd5f file
ds = ReactionLab.ModelData.Dataset.hdf5read(lc);

%% 
nUnit = length(ds.Targets);
dsName = ds.Title;
dsNew = generateDataset(dsName);
allID = {ds.OptimizationVariables.PrimeId};
allName = {ds.OptimizationVariables.Key};
LBall = [ds.OptimizationVariables.LowerBound]';
UBall = [ds.OptimizationVariables.UpperBound]';
Valall = [ds.OptimizationVariables.Value]';
dsNameAll = {ds.Targets.PrimeId};
dsBDall = reshape([ds.Targets.Bounds],2,nUnit);
dsBDall = dsBDall';
dsValall = [ds.Targets.Value]';

%% Set variables
sm = ds.SurrogateModels;
for i = 1:nUnit
   smID = {sm(i).Target.primeId};
   coefMatrix = sm(i).Coef;
   if isempty([sm(i).OptimizationVariables.varId])
      varID = {sm(i).OptimizationVariables.varPrimeId};
   else
      varID = {sm(i).OptimizationVariables.varId};
   end
   [~,~,idx] = intersect(varID,allID,'stable');
   varName = allName(idx);
   varLB = LBall(idx);
   varUB = UBall(idx);
   varVal = Valall(idx);
   H = [varLB, varUB];
   vars = generateVar(varName, H, varVal);
   Smodel = generateModel(coefMatrix,vars);
   dsUnit = generateDSunit(dsNameAll{i},Smodel,[dsBDall(i,1),dsBDall(i,2)],dsValall(i));
   dsNew.addDSunit(dsUnit);
end
   