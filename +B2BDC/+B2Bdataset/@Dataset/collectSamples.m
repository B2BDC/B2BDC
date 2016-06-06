function Xvals = collectSamples(obj,N)
%   XVALS = COLLECTSAMPLES(OBJ, N) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of 
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

opt = generateOpt('Display',false','ExtraLinFraction',0);
if obj.isConsistent(opt)
   xInit = obj.FeasiblePoint;
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end

% Check Dataset model type

Xvals = quadSamples(obj,xInit,N);


function Xvals = quadSamples(obj,xInit,N)
nTargets = obj.Length;
nX = numel(obj.VarNames);

Mgd = cell(2*nTargets,1);
idx = cell(2*nTargets,1);
for i=1:nTargets
   EUi = obj.DatasetUnits.Values(i);
   % Use of STABLE below ultimately preserves order in E.VarNames
   [~,~,VIDX] = intersect(EUi.SurrogateModel.VarNames, obj.VarNames, 'stable');
   idx{i} = VIDX;
   idx{i+nTargets} = VIDX;
   if isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.QModel')
      Mgd{i} = EUi.SurrogateModel.CoefMatrix;
      Mgd{i}(1,1) = Mgd{i}(1,1) - EUi.UpperBound;
      Mgd{i+nTargets} = -EUi.SurrogateModel.CoefMatrix;
      Mgd{i+nTargets}(1,1) = Mgd{i+nTargets}(1,1) + EUi.LowerBound;
   elseif isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.RQModel')
      Ni = EUi.SurrogateModel.Numerator;
      Di = EUi.SurrogateModel.Denominator;
      Mgd{i} = Ni - EUi.UpperBound*Di;
      Mgd{i+nTargets} = EUi.LowerBound*Di - Ni;
   else
      error('The model type is invalid')
   end
end

H = obj.Variables.calBound;

Xvals = zeros(N, nX);
Xvals(1,:) = xInit.';
for i=2:N
   V = randn(1, nX);
   zPDir = find(Xvals(i-1,:)' >= H*[0.01; 0.99]);
   zNDir = find(Xvals(i-1,:)' <= H*[0.99; 0.01]);
   V(zPDir) = -abs(V(zPDir));
   V(zNDir) = abs(V(zNDir));
   tmp = B2BDC.B2Bdataset.Dataset.q2sample(Mgd,idx,H,Xvals(i-1,:)',V');
   Xvals(i,:) = tmp';
end
