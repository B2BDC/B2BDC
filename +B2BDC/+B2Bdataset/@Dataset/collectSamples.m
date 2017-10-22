function Xvals = collectSamples(obj,N,x0,opt)
%   XVALS = COLLECTSAMPLES(OBJ, N, X0,OPT) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

if nargin < 4
   opt = generateOpt('Display',false);
end
if obj.FeasibleFlag
   opt.AddFitError = true;
end
if obj.isConsistent(opt)
   if nargin < 3
      xInit = obj.FeasiblePoint;
   else
      if isempty(x0)
         xInit = obj.FeasiblePoint;
      elseif obj.isFeasiblePoint(x0)
         if size(x0,1) == 1
            xInit = x0';
         else
            xInit = x0;
         end 
      else
         error('The input starting point is not feasible');
      end
   end
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end
Sopt = opt.SampleOption;
units = obj.DatasetUnits.Values;
abE = zeros(length(units),1);
if obj.FeasibleFlag
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end

% Check Dataset model type

Xvals = quadSamples(obj,xInit,N);


   function Xvals = quadSamples(obj,xInit,N)
      nTargets = obj.Length;
      nX = numel(obj.VarNames);
      vList = obj.Variables;
      A0 = vList.ExtraLinConstraint.A;
      if ~isempty(A0)
         nVar = vList.Length;
<<<<<<< Updated upstream
         tolerance = 1e-5;
         A0 = vList.ExtraLinConstraint.A;
         Aub = vList.ExtraLinConstraint.UB;
         Alb = vList.ExtraLinConstraint.LB;
         eqTest = (Aub-Alb)./sum(abs(A0),2);
         ieq = (eqTest <= tolerance);
         A0 = A0(~ieq,:);
         Aub = Aub(~ieq);
         Alb = Alb(~ieq);
=======
%          tolerance = 1e-5;
%          A0 = vList.ExtraLinConstraint.A;
         Aub = vList.ExtraLinConstraint.UB;
         Alb = vList.ExtraLinConstraint.LB;
%          eqTest = (Aub-Alb)./sum(abs(A0),2);
%          ieq = (eqTest <= tolerance);
%          A0 = A0(~ieq,:);
%          Aub = Aub(~ieq);
%          Alb = Alb(~ieq);
>>>>>>> Stashed changes
         nL = size(A0,1);
      else
         nL = 0;
      end
      
      
      Mgd = cell(2*nTargets+nL,1);
      idx = cell(2*nTargets+nL,1);
      for i=1:nTargets
         EUi = obj.DatasetUnits.Values(i);
         % Use of STABLE below ultimately preserves order in E.VarNames
         [~,~,VIDX] = intersect(EUi.SurrogateModel.VarNames, obj.VarNames, 'stable');
         idx{i} = VIDX;
         idx{i+nTargets} = VIDX;
         if isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.QModel')
            Mgd{i} = EUi.SurrogateModel.CoefMatrix;
            Mgd{i}(1,1) = Mgd{i}(1,1) - EUi.UpperBound -abE(i);
            Mgd{i+nTargets} = -EUi.SurrogateModel.CoefMatrix;
            Mgd{i+nTargets}(1,1) = Mgd{i+nTargets}(1,1) + EUi.LowerBound - abE(i);
         elseif isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            Ni = EUi.SurrogateModel.Numerator;
            Di = EUi.SurrogateModel.Denominator;
            Mgd{i} = Ni - (EUi.UpperBound+abE(i))*Di;
            Mgd{i+nTargets} = (EUi.LowerBound-abE(i))*Di - Ni;
         else
            error('The model type is invalid')
         end
      end
      
      for i = 1:nL
         ai = A0(i,:);
         ub = Aub(i);
         lb = Alb(i);
         tmpM = zeros(nVar+1);
         tmpM(2:end,2:end) = ai' * ai;
         tmpM(1) = ub*lb;
         tmpM(1,2:end) = -0.5*(ub+lb)*ai;
         tmpM(2:end,1) = -0.5*(ub+lb)*ai';
         Mgd{i+2*nTargets} = tmpM;
         idx{i+2*nTargets} = (1:nVar)';
      end
      
      H = obj.Variables.calBound;
      
      Xvals = zeros(N, nX);
      ic = 1;
%       h = waitbar((ic-1)/N,'Collecting samples in the feasible set...');
      nstep = Sopt.StepInterval;
      n2 = nstep*N;
      for i=1:n2
         V = randn(1, nX);
         zPDir = find(xInit >= H*[0.01; 0.99]);
         zNDir = find(xInit <= H*[0.99; 0.01]);
         V(zPDir) = -abs(V(zPDir));
         V(zNDir) = abs(V(zNDir));
         tmpx = B2BDC.B2Bdataset.Dataset.q2sample(Mgd,idx,H,xInit,V');
         if mod(i,nstep) == 0
            Xvals(ic,:) = tmpx';
            ic = ic+1;
%             if mod(ic,round(0.01*N)) == 1
%                waitbar((ic-1)/N,h);
%             end
         end
         xInit = tmpx;
      end
%       close(h);
%       Xvals = Xvals(randperm(N,N),:);
   end

end