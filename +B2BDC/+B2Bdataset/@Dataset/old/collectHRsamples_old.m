function xSample = collectHRsamples(obj,N,x0,opt)
%   XHR = COLLECTHRSAMPLES(OBJ, N, X0,OPT) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

%  Created: June 28, 2018     Wenyu Li

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
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end
tolX = 1e-3;
tolY = 1e-4;
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
% ddy = 0.5*tolY*(UB-LB);
% UB = UB;
% LB = LB;
nStep = opt.SampleOption.StepInterval;
d = randn(nVar+nMD+nPD,nStep*N);
d = normc(d);
% ts = (1-0.5*tolerance)*rand(nStep*N,1)+0.5*tolerance;
ts = rand(nStep*N,1);
Ad = Ax*d;
xx = xInit;
xHR = zeros(N,nVar+nMD+nPD);
nS = 1;
if opt.Display
   h = waitbar(0,'Collecting hit and run samples...');
   for i = 1:nStep*N
      tt = calculateIntersection(d(:,i),xx);
      xx = xx + tt*d(:,i);
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
      if mod(nS,round(0.05*N)) == 0
         waitbar(nS/N,h);
      end
   end
else
   yQOI = [];
   for i = 1:nStep*N
      tt = calculateIntersection(d(:,i),xx);
      xx = xx + tt*d(:,i);
      if mod(i,nStep) == 0
         if ~obj.isFeasiblePoint(xx')
            keyboard;
         end
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


   function tt = calculateIntersection(d,xx)
      tb = bx-Ax*xx;
      tv = Ad(:,i);
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
         vv = d(idx{ii});
         aa(ii) = vv'*Coef(2:end,2:end)*vv;
         bb(ii) = 2*(vv'*Coef(2:end,2:end)*xx(idx{ii})+Coef(1,2:end)*vv);
         cc = [1;xx(idx{ii})]'*Coef*[1;xx(idx{ii})];
         cUB(ii) = cc-UB(ii);
         cLB(ii) = cc-LB(ii);
      end
      f = @QuickEval;
      delta_UB = bb.^2-4*aa.*cUB;
      delta_LB = bb.^2-4*aa.*cLB;
      flag_UB = find(delta_UB < 0);
      flag_LB = find(delta_LB < 0);
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
      if nP>1
         tQPm = 0.5*(tQP(1:end-1)+tQP(2:end));
      else
         tQPm = 0.5*tLP;
      end
      tQN = sort(tAll(tAll<=0),'descend');
      tQN(tQN<=tLN) = [];
      nN = length(tQN);
      tQN(end+1,1) = tLN;
      if nN>1
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
               if tQ(ti,1) > tQ(ti,tj) && tQ(ti,1) < tQ(ti,4)
                  idp = find(tQP == tQ(ti,1),1);
               else
                  idp = find(tQP == tQ(ti,4),1);
               end
            elseif tj == 2
               if tQ(ti,4) > tQ(ti,tj)
                  idp = find(tQP == tQ(ti,4),1);
               else
                  break
               end
            end
         else
            if tj == 1
               if tQ(ti,3) > tQ(ti,tj) && tQ(ti,3) < tQ(ti,2)
                  idp = find(tQP == tQ(ti,3),1);
               else
                  idp = find(tQP == tQ(ti,2),1);
               end
            elseif tj == 4
               if tQ(ti,2) > tQ(ti,tj)
                  idp = find(tQP == tQ(ti,2),1);
               else
                  break
               end
            end
         end
         for jj = idp:nP-2
            yy = f(tQPm(jj));
            if all(yy<0)
               idp = jj+1;
               tp = [tp; tQP(idp-1:idp)'];
               break
            end
         end
         yy = f(tQPm(nP-1));
         if all(yy<0)
            tp = [tp; tQP(nP-1:nP)'];
         end
         break
      end
      if nP > 0
         yy = f(tQPm(nP));
         if all(yy<0)
            tp = [tp; tQP(end-1:end)'];
         end
      end
      % negative t
      tn = [max(tQN) 0];
      idn = 1;
      while idn < nN
         [ti,tj] = find(tQ==tn(end,1));
         if aa(ti) > 0
            if tj == 4
               if tQ(ti,2) < tQ(ti,tj) && tQ(ti,2) > tQ(ti,3)
                  idn = find(tQN == tQ(ti,2),1);
               else
                  idn = find(tQN == tQ(ti,3),1);
               end
            elseif tj == 1
               if tQ(ti,3) < tQ(ti,tj)
                  idn = find(tQN == tQ(ti,3),1);
               else
                  break
               end
            end
         else
            if tj == 2
               if tQ(ti,4) < tQ(ti,tj) && tQ(ti,4) > tQ(ti,1)
                  idn = find(tQN == tQ(ti,4),1);
               else
                  idn = find(tQP == tQ(ti,1),1);
               end
            elseif tj == 3
               if tQ(ti,1) < tQ(ti,tj)
                  idn = find(tQP == tQ(ti,1),1);
               else
                  break
               end
            end
         end
         for jj = idn:nN-2
            yy = f(tQNm(jj));
            if all(yy<0)
               idn = jj+1;
               tn = [tn; tQN([idn,idn-1])'];
               break
            end
         end
         yy = f(tQNm(nN-1));
         if all(yy<0)
            tn = [tn; tQN([nN,nN-1])'];
         end
         break
      end
      if nN > 0
         yy = f(tQNm(nN));
         if all(yy<0)
            tn = [tn; tQN([nN+1,nN])'];
         end
      end
      tcand = [tp; tn];
      dt = diff(tcand,[],2);
      tsum = [0;cumsum(dt)/sum(dt)];
      it = find(tsum>ts(i),1);
      frac = (ts(i)-tsum(it-1))/(tsum(it)-tsum(it-1));
      % not too close to the boundary
      if frac < 0.5*tolX
         frac = 0.5*tolX;
      elseif frac > 1-0.5*tolX
         frac = 1-0.5*tolX;
      end
      tt = tcand(it-1,1)+dt(it-1)*frac;
        
         function y = QuickEval(t)
            y = zeros(2*nUnit,1);
            y(1:nUnit) = aa*t^2+bb*t+cUB;
            y(nUnit+1:end) = -(aa*t^2+bb*t+cLB);
         end
   end
end