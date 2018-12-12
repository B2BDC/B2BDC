function [xSample,flag,CC] = collectHRsamples_BCW(obj,N,x0,opt,nB)
%   XHR = COLLECTHRSAMPLES(OBJ, N, X0,OPT) returns a n-by-nVariable matrix of
%   feasible points of the dataset OBJ, where nVariable is the number of
%   variables of the B2BDC.B2Bdataset.Dataset object and N is the number of
%   feasible points to be returned.

%  Created: June 28, 2018     Wenyu Li
%  Modified: Oct 15, 2018     Wenyu Li

flag = 0;
CC = 0;
tolX = 0;
% tolY = 1e-4;
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
sB = ceil(nV/nB);
ts = rand(nStep*N*sB,1);
aFlag = cell(sB,1);
tmpID = 1:nUnit;
d = cell(sB,1);
for i1 = 1:sB
   if i1 == sB && mod(nV,nB) ~= 0
      tID = nV-mod(nV,nB)+1:nV;
      dd = randn(mod(nV,nB),nStep*N);
   else
      tID = (1:nB)+nB*(i1-1);
      dd = randn(nB,nStep*N);
   end
   d{i1} = normc(dd);
   tmpFlag = true(nUnit,1);
   for i2 = 1:nUnit
      [~,iid] = intersect(idx{i2},tID);
      if isempty(iid)
         tmpFlag(i2) = false;
      end
   end
   aFlag{i1} = tmpID(tmpFlag);
end
xx = xInit;
xHR = zeros(N,nV);
nS = 1;
if opt.Display
   h = waitbar(0,'Collecting hit and run samples...');
   for i = 1:nStep*N
      for j = 1:sB
         if j == sB && mod(nV,nB) ~= 0
            tID = nV-mod(nV,nB)+1:nV;
         else
            tID = (1:nB)+(j-1)*nB;
         end
         [tt,tcc] = calculateIntersection(xx);
         CC = CC+tcc;
         xx(tID) = xx(tID) + tt*d{j}(:,i);
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
      for j = 1:sB
         if j == sB && mod(nV,nB) ~= 0
            tID = nV-mod(nV,nB)+1:nV;
         else
            tID = (1:nB)+(j-1)*nB;
         end
         [tt,tcc] = calculateIntersection(xx);
         CC = CC + tcc;
         xx(tID) = xx(tID) + tt*d{j}(:,i);
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
CC = CC/nStep/N/sB;
% if sum(obj.isFeasiblePoint(xHR)) < N
%    disp('Errors');
% end
xSample.x = xHR;
xSample.dimension = [nVar nMD nPD];


   function [tt,count] = calculateIntersection(xx)
      D = zeros(nV,1);
      D(tID) = d{j}(:,i);
      tb = bx-Ax*xx;
      tv = Ax*D;
      Pos = tv>0;
      Neg = tv<0;
      tLP = min(tb(Pos)./tv(Pos));
      tLN = max(tb(Neg)./tv(Neg));
      aa = inf(nUnit,1);
      bb = zeros(nUnit,1);
      cUB = ones(nUnit,1);
      cLB = ones(nUnit,1);
      for ii = aFlag{j}
         Coef = Qx{ii};
         vv = D(idx{ii});
         aa(ii) = vv'*Coef(2:end,2:end)*vv;
         bb(ii) = 2*(vv'*Coef(2:end,2:end)*xx(idx{ii})+Coef(1,2:end)*vv);
         cc = [1;xx(idx{ii})]'*Coef*[1;xx(idx{ii})];
         cUB(ii) = cc-UB(ii);
         cLB(ii) = cc-LB(ii);
      end
      f = @(t) QuickEval(t,aFlag{j});
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
      count = 0;
      % positive t
      tp = [0 0];
      idp = 1;
      if nP == 0
         tp(end,2) = tQP;
      else
         while idp <= nP
            % search for infeasible right end
            while idp <= nP && f(tQPm(idp))<0
               idp = idp+1;
               count = count+1;
            end
            count = count+1;
            if idp > nP
               tp(end,2) = tQP(end);
               break
            else
               tp(end,2) = tQP(idp);
            end
            % find next candidate feasible left end point
            [ti,tj] = find(tQ==tp(end,2));
            if aa(ti) > 0
               if tj == 3
                  if tQ(ti,1) > tQ(ti,tj) && tQ(ti,1) < tQ(ti,4) && ~isempty(find(tQP == tQ(ti,1),1))
                     idp = find(tQP == tQ(ti,1),1);
                  elseif ~isempty(find(tQP == tQ(ti,4),1))
                     idp = find(tQP == tQ(ti,4),1);
                  else
                     break;
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
            % find a feasible left end point
            id0 = idp;
            for jj = idp:nP-1
               yy = f(tQPm(jj));
               count = count+1;
               if yy<0
                  tp(end+1,1) = tQP(jj);
                  idp = jj+1;
                  break
               end
            end
            if idp > id0
               continue;
            elseif f(tQPm(end))<0
               tp(end+1,:) = tQP(end-1:end);
               count = count+1;
               break
            else
               break
            end
         end
      end
      
      % negative t
      tn = [0 0];
      idn = 1;
      if nN == 0
         tn(end,1) = tQN;
      else
         while idn <= nN
            % search for infeasible left end
            while idn <= nN && f(tQNm(idn))<0
               idn = idn+1;
               count = count+1;
            end
            count = count+1;
            if idn > nN
               tn(end,1) = tQN(end);
               break
            else
               tn(end,1) = tQN(idn);
            end
            % find next candidate feasible right end point
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
            % find a feasible right end point
            id0 = idn;
            for jj = idn:nN-1
               yy = f(tQNm(jj));
               count = count+1;
               if yy<0
                  tn(end+1,2) = tQN(jj);
                  idn = jj+1;
                  break
               end
            end
            if idn > id0
               continue
            elseif f(tQNm(end))<0
               tn(end+1,:) = flipud(tQN(end-1:end));
               count = count+1;
               break
            else
               break
            end
         end
      end
      tn(1,2) = tp(1,2);
      tcand = [tp(2:end,:); tn];
      dt = diff(tcand,[],2);
      tsum = [0;cumsum(dt)/sum(dt)];
      it = find(tsum>ts((i-1)*(sB)+j),1);
      if ~isempty(it)
         frac = (ts((i-1)*(sB)+j)-tsum(it-1))/(tsum(it)-tsum(it-1));
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
      
      function y = QuickEval(t,flag)
         nn = length(flag);
         y = repmat(aa(flag)*t^2+bb(flag)*t,2,1)+[cUB(flag); cLB(flag)];
         y(nn+1:end) = -y(nn+1:end);
         y = max(y);
      end
   end
end