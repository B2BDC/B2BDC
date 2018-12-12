function setModelDiscrepancy(obj,QOIindex,fMD,Range)
% SETMODELDISCREPANCY(OBJ,QOIINDEX,SV,FMD,SVRANGE) sets the model discrepancy
% correction to the dataset object.

% obj - B2BDC.B2Bdataset.Dataset object
% QOIindex - A cell array of groups of QOIs to add model discrepancy correction
% fMD - A cell array of function handles that generates basis functions for model discrepancy correction
% Range -  A numerical array of maginitude constraint on correction coefficients

%     Created: Sep 5, 2018     Wenyu Li

obj.clearModelDiscrepancy;
if isempty(QOIindex)
   QOIindex{1} = 1:obj.Length;
end
if ~iscell(QOIindex)
   QOIindex = {QOIindex};
end
if ~iscell(fMD)
   fMD = {fMD};
end
nGroup = length(QOIindex);
if nargin < 4 || isempty(Range)
   Range = repmat(10,nGroup,1);
elseif isscalar(Range)
   Range = repmat(Range,nGroup,1);
elseif length(Range) ~= nGroup || any(Range <= 0)
   error('The input range of scenario parameter is invalid')
end
if length(fMD) ~= nGroup
   error('The input basis function has a wrong dimension')
end
nQOI = obj.Length;
nmd = zeros(1,nGroup);
sv = obj.getScenarioParameter;
GroupIndex = zeros(nQOI,1);
v0 = B2BDC.B2Bvariables.VariableList;
Basis = cell(nQOI,1);
for i = 1:nGroup
   qid = QOIindex{i};
   if any(qid>obj.Length)
      error('Input QOI index exceeds dataset length');
   end
   GroupIndex(qid) = i;
   basis = fMD{i}(sv(qid,:));
   nmd(i) = size(basis,2);
   vname = cell(nmd(i),1);
   for j = 1:nmd(i)
      vname{j} = ['Correction coef. ' num2str(j) ' (Group: ' num2str(i) ')'];
   end
   vGroup = generateVar(vname,repmat([-1 1]*Range(i),nmd(i),1));
   v0 = v0.addList(vGroup);
   for j = 1:length(qid)
      Basis{qid(j)} = basis(j,:);
   end
end
rr.Variables = v0;
rr.GroupIndex = GroupIndex;
rr.CorrectionDimension = nmd;
rr.Basis = Basis;
rr.FeasiblePoint = [];
rr.BasisFunction = fMD;
obj.ModelDiscrepancy = rr;
obj.ModelDiscrepancyFlag = true;

