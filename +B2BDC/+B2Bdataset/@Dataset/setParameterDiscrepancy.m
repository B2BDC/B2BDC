function setParameterDiscrepancy(obj,Varindex,fPD,Hv,Range)
% SETMODELDISCREPANCY(OBJ,QOIINDEX,SV,FMD,SVRANGE) sets the model discrepancy
% correction to the dataset object.

% obj - B2BDC.B2Bdataset.Dataset object
% Varindex - A numerical array specifies variables to add parameter discrepancy correction
% fPD - A cell array of function handles that generates basis functions for parameter discrepancy correction
% Hv - A numerical array specifies the range for the corrected model parameters
% Range -  A numerical array of maginitude constraint on correction coefficients

%     Created: Sep 7, 2018     Wenyu Li

obj.clearParameterDiscrepancy;
[Varindex,II] = sort(Varindex,'ascend');
nVar = obj.Variables.Length;
if any(Varindex>nVar)
   error('Input variable index exceeds dataset variable length');
end
nCorrect = length(Varindex);
if nargin < 4
   Hv = [];
end
if isempty(Hv)
   Hv = repmat([-inf inf],nVar,1);
elseif size(Hv,1) == nCorrect
   hv = repmat([-inf inf],nVar,1);
   hv(Varindex,:) = Hv;
   Hv = hv;
end
if size(Hv,1) ~= nVar
   error('Wrong dimension of input model parameter range');
end
if ~iscell(fPD)
   fPD = {fPD};
end
fPD = fPD(II);
if nargin < 5 || isempty(Range)
   Range = repmat(10,nCorrect,1);
elseif isscalar(Range)
   Range = repmat(Range,nCorrect,1);
elseif length(Range) ~= nCorrect || any(Range <= 0)
   error('The input range of scenario parameter is invalid')
end
if length(fPD) ~= nCorrect
   error('The input basis function has a wrong dimension')
end
npd = zeros(1,nVar);
sv = obj.getScenarioParameter;
v0 = B2BDC.B2Bvariables.VariableList;
for i = 1:nCorrect
   vid = Varindex(i);
   npd(vid) = length(fPD{i}(sv(1,:)));
   vname = cell(npd(vid),1);
   for j = 1:npd(vid)
      vname{j} = ['Correction coef. ' num2str(j) ' (Variable: ' num2str(vid) ')'];
   end
   vcorrect = generateVar(vname,repmat([-1 1]*Range(i),npd(vid),1));
   v0 = v0.addList(vcorrect);
end
rr.Variables = v0;
rr.CorrectionDimension = npd;
nQOI = obj.Length;
vNames = obj.VarNames;
Basis = cell(nQOI,1);
for i = 1:nQOI
   tmpModel = obj.DatasetUnits.Values(i).SurrogateModel;
   tmpName = tmpModel.VarNames;
   basis = cell(length(tmpName),1);
   [~,~,id1] = intersect(tmpName,vNames,'stable');
   for j = 1:length(id1)
      id2 = find(Varindex == id1(j));
      if ~isempty(id2)
         basis{j} = fPD{id2}(sv(i,:));
      end
   end
   Basis{i} = basis;
end
rr.Basis = Basis;
rr.VariableRange = Hv;
rr.FeasiblePoint = [];
rr.BasisFunction = fPD;
obj.ParameterDiscrepancy = rr;
obj.ParameterDiscrepancyFlag = true;