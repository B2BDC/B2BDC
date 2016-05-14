function dc = datasetToDCContainer(obj)
% This function create a DClab.DCContainer object from a
% B2BDC.B2Bdataset.Dataset object.

%  Created: Nov 26, 2015     Wenyu Li

dc = DClab.DCContainer;
dc.title = obj.Name;
dc.parameterKey = obj.VarNames;
nVar = obj.Variables.Length;
tep_varRange = cell(nVar,1);
varRange = obj.Variables.calBound;
for i = 1:nVar
   tep_varRange{i} = varRange(i,:);
end
dc.parameterRange = tep_varRange;
tep_parID = cell(nVar,1);
for i = 1:nVar
   tep_parID{i} = ['v' num2str(i)];
end
dc.parameterID = tep_parID;
tep_Initial = cell(nVar,1);
tep_ob = [obj.Variables.Values.NominalValue];
for i = 1:nVar
   tep_Initial{i} = tep_ob(i);
end
dc.paramInitialValue = tep_Initial;
nUnit = obj.Length;
dc.modelKey = {obj.DatasetUnits.Values.Name}';
tep_coef = cell(nUnit,1);
for i = 1:nUnit
   tep_coef{i} = obj.DatasetUnits.Values(i).SurrogateModel.NormalizedCoefMatrix;
end
dc.modelCoeffs = tep_coef;

tep_modelParamID = cell(nUnit,1);
for i = 1:nUnit
   [~,~,id] = intersect(obj.DatasetUnits.Values(i).SurrogateModel.VarNames,obj.VarNames,'stable');
   tep_modelParamID{i} = tep_parID(id);
end
dc.modelParamIDs = tep_modelParamID;

dc.targetKey = {obj.DatasetUnits.Values.Name}';
expBound = obj.calBound;
tep_uncRange = cell(nUnit,1);
for i = 1:nUnit
   tep_uncRange{i} = expBound(i,:);
end
dc.uncRange = tep_uncRange;
dc.targetValue = {obj.DatasetUnits.Values.ObservedValue}';
