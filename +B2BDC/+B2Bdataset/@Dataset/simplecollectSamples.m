function [Xvals] = simplecollectSamples(E,N)
% Returns N feasible points of the dataset. The returned Xvals is a
% nSample-by-nVariable matrix.

opt = generateOpt('Display',false','ExtraLinFraction',0);
if E.isConsistent(opt)
   xInit = E.FeasiblePoint;
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end

% Check Dataset model type

Xvals = quadSamples(E,xInit,N);


function Xvals = quadSamples(E,xInit,N)
nTargets = E.Length;
nX = numel(E.VarNames);

Mgd = cell(2*nTargets,1);
idx = cell(2*nTargets,1);
for i=1:nTargets
   EUi = E.DatasetUnits.Values(i);
   % Use of STABLE below ultimately preserves order in E.VarNames
   [~,~,VIDX] = intersect(EUi.SurrogateModel.VarNames, E.VarNames, 'stable');
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

H = E.Variables.calBound;

Xvals = zeros(N, nX);
Xvals(1,:) = xInit.';
for i=2:N
   V = randn(1, nX);
   zPDir = find(Xvals(i-1,:)' >= H*[0.01; 0.99]);
   zNDir = find(Xvals(i-1,:)' <= H*[0.99; 0.01]);
   V(zPDir) = -abs(V(zPDir));
   V(zNDir) = abs(V(zNDir));
   tmp = q2sample(Mgd,idx,H,Xvals(i-1,:)',V');
   Xvals(i,:) = tmp';
end
