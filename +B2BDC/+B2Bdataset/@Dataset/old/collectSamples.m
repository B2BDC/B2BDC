function xSample = collectSamples(obj,N,x0,opt)
%   XVALS = COLLECTSAMPLES(OBJ, N, X0,OPT) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

if nargin < 4
   opt = generateOpt('Display',false,'ConsistencyMeasure','absolute');
end
if nargin < 3
   x0 = [];
end
if obj.FeasibleFlag
   opt.AddFitError = true;
end
if ~obj.isConsistent(opt)
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

xSample = quadSamples(obj,x0,N);


   function xSample = quadSamples(obj,x0,N)
      if obj.ModelDiscrepancyFlag
         nMD = obj.ModelDiscrepancy.Variables.Length;
      else
         nMD = 0;
      end
      if obj.ParameterDiscrepancyFlag
         nPD = obj.ParameterDiscrepancy.Variables.Length;
      else
         nPD = 0;
      end
      if isempty(x0)
         xInit = obj.FeasiblePoint;
         if nMD > 0
            xInit = [xInit; obj.ModelDiscrepancy.FeasiblePoint];
         end
         if nPD > 0
            xInit = [xInit; obj.ParameterDiscrepancy.FeasiblePoint];
         end
      else
         if size(x0,1) == 1
            x0 = x0';
         end
         [~,xInit] = obj.isFeasiblePoint(x0');
         if isempty(xInit)
            error('The input starting point is not feasible');
         else
            xInit = xInit';
         end
      end
      nTargets = obj.Length;
      vList = obj.Variables;
      nVar = vList.Length;
      nX = nVar+nMD+nPD;
      [idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion;
      A0 = vList.ExtraLinConstraint.A;
      if ~isempty(A0)
         A0 = [A0 zeros(size(A0,1), nMD+nPD)];
         Aub = vList.ExtraLinConstraint.UB;
         Alb = vList.ExtraLinConstraint.LB;
      else
         Aub = [];
         Alb = [];
      end
      A0 = [A0; APD(1:2:end,:)];
      Aub = [Aub; bPD(1:2:end)];
      Alb = [Alb; -bPD(2:2:end)];
      nL = size(A0,1);
 
      Mgd = cell(2*nTargets+nL,1);
      idx = cell(2*nTargets+nL,1);
      for i=1:nTargets
         EUi = obj.DatasetUnits.Values(i);
         idx{i} = idall{i};
         idx{i+nTargets} = idall{i};
         if isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.QModel')
            Mgd{i} = Qall{i};
            Mgd{i}(1,1) = Mgd{i}(1,1) - EUi.UpperBound -abE(i);
            Mgd{i+nTargets} = -Qall{i};
            Mgd{i+nTargets}(1,1) = Mgd{i+nTargets}(1,1) + EUi.LowerBound - abE(i);
         elseif isa(EUi.SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            Ni = Nall{i};
            Di = Dall{i};
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
         tmpM = zeros(nVar+nMD+nPD+1);
         tmpM(2:end,2:end) = ai' * ai;
         tmpM(1) = ub*lb;
         tmpM(1,2:end) = -0.5*(ub+lb)*ai;
         tmpM(2:end,1) = -0.5*(ub+lb)*ai';
         Mgd{i+2*nTargets} = tmpM;
         idx{i+2*nTargets} = (1:nVar+nMD+nPD)';
      end
      
      H = obj.Variables.calBound;
      if nMD > 0
         H = [H; obj.ModelDiscrepancy.Variables.calBound];
      end
      if nPD > 0
         H = [H; obj.ParameterDiscrepancy.Variables.calBound];
      end
      
      Xvals = zeros(N, nX);
      ic = 1;
      if opt.Display
         h = waitbar((ic-1)/N,'Collecting samples in the feasible set...');
      end
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
            if opt.Display
               if mod(ic,round(0.02*N)) == 1
                  waitbar((ic-1)/N,h);
               end
            end
         end
         xInit = tmpx;
      end
      if opt.Display
         close(h);
      end
      xSample.x = Xvals;
      xSample.dimension = [nVar nMD nPD];
%       Xvals = Xvals(randperm(N,N),:);
   end

end