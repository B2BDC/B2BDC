function [xSample,flag] = collectHRsamples_CW(obj,N,x0,opt)
%   XHR = COLLECTHRSAMPLES(OBJ, N, X0,OPT) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

%  Created: June 28, 2018     Wenyu Li
%  Modified: Oct 15, 2018     Wenyu Li

flag = 0;
if ~quadratictest(obj)
   error('HR sampling is now only available for quadratic models')
end
if nargin < 4
   opt = generateOpt('Display',false,'ConsistencyMeasure','absolute');
end
if nargin < 3
   x0 = [];
end
if obj.FeasibleFlag
   opt.AddFitError = true;
end
if obj.isConsistent(opt)
   if obj.ModelDiscrepancyFlag
      nMD = length(obj.ModelDiscrepancy.FeasiblePoint);
   else
      nMD = 0;
   end
   if obj.ParameterDiscrepancyFlag
      nPD = length(obj.ParameterDiscrepancy.FeasiblePoint);
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
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end
tolX = 0;
% tolY = 1e-4;
[idx,Qx,~,~,APD,bPD] = obj.getQ_RQ_expansion;
bounds = obj.calBound;
LB = bounds(:,1);
UB = bounds(:,2);
sVar = obj.Variables;
bounds = sVar.calBound;
if nMD > 0
   bounds = [bounds; obj.ModelDiscrepancy.Variables.calBound];
end
if nPD > 0
   bounds = [bounds; obj.ParameterDiscrepancy.Variables.calBound];
end
nVar = sVar.Length;
Ax = [eye(nVar+nMD+nPD); -eye(nVar+nMD+nPD)];
bx = [bounds(:,2); -bounds(:,1)];
if ~isempty(sVar.ExtraLinConstraint.A)
   AA = [sVar.ExtraLinConstraint.A zeros(size(sVar.ExtraLinConstraint.A,1),nMD+nPD)];
   Ax = [Ax; AA; -AA];
   bx = [bx; sVar.ExtraLinConstraint.UB; -sVar.ExtraLinConstraint.LB];
end
Ax = [Ax; APD];
bx = [bx; bPD];
units = obj.DatasetUnits.Values;
nUnit = length(units);
if obj.FeasibleFlag
   for i = 1:nUnit
      if ~isempty(units(i).SurrogateModel.ErrorStats)
         LB(i) = LB(i)-units(i).SurrogateModel.ErrorStats.absMax;
         UB(i) = UB(i)+units(i).SurrogateModel.ErrorStats.absMax;
      end
   end
end
idQ = [];
AQ = [];
nV = nVar+nMD+nPD;
for i = 1:nUnit
   tCoef = Qx{i}(2:end,2:end);
   if ~any(tCoef(:))
      tA = zeros(1,nV);
      tA(idx{i}) = 2*tCoef(1,2:end);
      AQ = [AQ; tA];
      idQ = [idQ;i];
   end
end
Ax = [Ax;AQ;-AQ];
bx = [bx;UB(idQ);-LB(idQ)];
idx(idQ) = [];
Qx(idQ) = [];
UB(idQ) = [];
LB(idQ) = [];
nUnit = length(UB);
% ddy = 0.5*tolY*(UB-LB);
% UB = UB;
% LB = LB;
nStep = opt.SampleOption.StepInterval;
ts = rand(nStep*N*nV,1);
xx = xInit;
xHR = zeros(N,nV);
nS = 1;
if opt.Display
   h = waitbar(0,'Collecting hit and run samples...');
   for i = 1:nStep*N
      for j = 1:nV
         tt = calculateIntersection(xx);
         xx(j) = xx(j) + tt;
      end
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
      if mod(nS,round(0.05*N)) == 0
         waitbar(nS/N,h);
      end
   end
else
%    yQOI = [];
   for i = 1:nStep*N
      for j = 1:nV
         tt = calculateIntersection(xx);
         xx(j) = xx(j) + tt;
      end
      if mod(i,nStep) == 0
%          if ~obj.isFeasiblePoint(xx')
%             keyboard;
%          end
         xHR(nS,:) = xx';
         nS = nS+1;
      end
   end
end
% if sum(obj.isFeasiblePoint(xHR)) < N
%    disp('Errors');
% end
xSample.x = xHR;
xSample.dimension = [nVar nMD nPD];


   function tt = calculateIntersection(xx)
      tb = bx-Ax*xx;
      tv = Ax(:,j);
      Pos = tv>0;
      Neg = tv<0;
      tLP = min(tb(Pos)./tv(Pos));
      tLN = max(tb(Neg)./tv(Neg));
      aa = zeros(nUnit,1);
      bb = zeros(nUnit,1);
      cUB = zeros(nUnit,1);
      cLB = zeros(nUnit,1);
      for ii = 1:nUnit
         Coef = Qx{ii};
         iid = find(idx{ii}==j,1);
         if ~isempty(iid)
            aa(ii) = Coef(iid+1,iid+1);
            bb(ii) = 2*(Coef(iid+1,2:end)*xx(idx{ii})+Coef(1,iid+1))-Coef(iid+1,iid+1)*xx(iid);
            cc = xx(iid)^2*Coef(iid+1,iid+1)+Coef(1,1)+2*Coef(1,iid+1)*xx(iid);
            cUB(ii) = cc-UB(ii);
            cLB(ii) = cc-LB(ii);
         else
            aa(ii) = inf;
            bb(ii) = 0;
            cUB(ii) = 1;
            cLB(ii) = 1;
         end
      end
      f = @QuickEval;
      delta_UB = bb.^2-4*aa.*cUB;
      delta_LB = bb.^2-4*aa.*cLB;
      flag_UB = find(delta_UB <= 0);
      flag_LB = find(delta_LB <= 0);
      aa = repmat(aa,1,4);
      bb = repmat(bb,1,4);
      cc = repmat([-1,1],nUnit,2).*sqrt([delta_UB delta_UB delta_LB delta_LB]);
      tQ = 0.5*(-bb + cc)./aa;
      tQ(flag_UB,1:2) = repmat([-inf inf],length(flag_UB),1);
      tQ(flag_LB,3:4) = repmat([-inf inf],length(flag_LB),1);
      tQ = real(tQ);
      aa = aa(:,1);
      bb = bb(:,1);
      tAll = tQ(:);
      tQP = sort(tAll(tAll>0),'ascend');
      tQP(tQP>=tLP) = [];
      nP = length(tQP);
      tQP(end+1,1) = tLP;
      if nP>0
         tQPm = 0.5*(tQP(1:end-1)+tQP(2:end));
      else
         tQPm = 0.5*tLP;
      end
      tQN = sort(tAll(tAll<=0),'descend');
      tQN(tQN<=tLN) = [];
      nN = length(tQN);
      tQN(end+1,1) = tLN;
      if nN>0
         tQNm = 0.5*(tQN(1:end-1)+tQN(2:end));
      else
         tQNm = 0.5*tLN;
      end
      % positive t
      tp = [0 min(tQP)];
      idp = 1;
      while idp < nP
         [ti,tj] = find(tQ==tp(end,2));
         if aa(ti) > 0
            if tj == 3
               if tQ(ti,1) > tQ(ti,tj) && tQ(ti,1) < tQ(ti,4) && ~isempty(find(tQP == tQ(ti,1),1))
                  idp = find(tQP == tQ(ti,1),1);
               elseif ~isempty(find(tQP == tQ(ti,4),1))
                  idp = find(tQP == tQ(ti,4),1);
               else
                  break
               end
            else
               break
            end
         else
            if tj == 1
               if tQ(ti,3) > tQ(ti,tj) && tQ(ti,3) < tQ(ti,2) && ~isempty(find(tQP == tQ(ti,3),1))
                  idp = find(tQP == tQ(ti,3),1);
               elseif ~isempty(find(tQP == tQ(ti,2),1))
                  idp = find(tQP == tQ(ti,2),1);
               else
                  break;
               end
            else
               break;
            end
         end
         yy = [];
         for jj = idp:nP-1
            yy = f(tQPm(jj));
            if all(yy<0)
               idp = jj+1;
               tp = [tp; tQP(jj:jj+1)'];
               break
            end
         end
         if isempty(yy) || any(yy>0)
            break
         end
      end
      yy = f(tQPm(end));
      if nP > 0 && all(yy<0)
         tp = [tp; tQP(end-1:end)'];
      end
      % negative t
      tn = [max(tQN) 0];
      idn = 1;
      while idn < nN
         [ti,tj] = find(tQ==tn(end,1));
         if aa(ti) > 0
            if tj == 4
               if tQ(ti,2) < tQ(ti,tj) && tQ(ti,2) > tQ(ti,3) && ~isempty(find(tQN == tQ(ti,2),1))
                  idn = find(tQN == tQ(ti,2),1);
               elseif ~isempty(find(tQN == tQ(ti,3),1))
                  idn = find(tQN == tQ(ti,3),1);
               else
                  break;
               end
            else
               break;
            end
         else
            if tj == 2
               if tQ(ti,4) < tQ(ti,tj) && tQ(ti,4) > tQ(ti,1) && ~isempty(find(tQN == tQ(ti,4),1))
                  idn = find(tQN == tQ(ti,4),1);
               elseif ~isempty(find(tQP == tQ(ti,1),1))
                  idn = find(tQP == tQ(ti,1),1);
               else
                  break
               end
            else
               break
            end
         end
         yy = [];
         for jj = idn:nN-1
            yy = f(tQNm(jj));
            if all(yy<0)
               idn = jj+1;
               tn = [tn; tQN([jj+1,jj])'];
               break
            end
         end
         if isempty(yy) || any(yy>0)
            break;
         end
      end
      yy = f(tQNm(end));
      if nN > 0 && all(yy<0)
         tn = [tn; tQN([end,end-1])'];
      end
      tcand = [tp; tn];
      tn(1,2) = tp(1,2);
      tcand = [tp(2:end,:); tn];
      dt = diff(tcand,[],2);
      tsum = [0;cumsum(dt)/sum(dt)];
      it = find(tsum>ts((i-1)*(nV)+j),1);
      if ~isempty(it)
         frac = (ts((i-1)*(nV)+j)-tsum(it-1))/(tsum(it)-tsum(it-1));
         % not too close to the boundary
         if frac < 0.5*tolX
            frac = 0.5*tolX;
         elseif frac > 1-0.5*tolX
            frac = 1-0.5*tolX;
         end
      else
         frac = 0;
         flag = flag+1;
      end
      tt = tcand(it-1,1)+dt(it-1)*frac;
        
         function y = QuickEval(t)
            y = zeros(2*nUnit,1);
            y(1:nUnit) = aa*t^2+bb*t+cUB;
            y(nUnit+1:end) = -(aa*t^2+bb*t+cLB);
         end
   end
end