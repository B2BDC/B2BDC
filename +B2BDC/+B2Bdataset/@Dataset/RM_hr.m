function [xVal,eff,opt] = RM_hr(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 3, 2016     Wenyu Li

% tic;
th = opt.SampleOption.RejectionTol;
p1 = 1-opt.SampleOption.PCTruncation;
xVal = [];
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nCut = sopt.ExtraCut;
if isempty(Dinfo)
   Dinfo.D = [];
   Dinfo.uq = [];
end
if (p1>0 && isempty(sopt.PCinfo)) ||...
      ((strcmp(s1,'Sample') || strcmp(s1,'Truncation')) && nCut >0)
   if isempty(xS)
      nS = min(10^7,10*nVar^3);
      xS = obj.collectSamples(nS,[],opt);
   else
      nS = size(xS,1);
   end
end
if p1 > 0
   if ~isempty(xS)
      xAve = mean(xS);
      xC = xS - repmat(xAve,nS,1);
      [V,DD] = eig((xC' * xC)/nS);
      Ddiag = diag(DD);
      [dd,id] = sort(Ddiag,'descend');
      dd = dd/sum(dd);
      vv = V(:,id);
      pcinfo.direction = vv;
      pcinfo.variance = dd;
      opt.SampleOption.PCinfo = pcinfo;
   else
      xAve = sopt.PCinfo.mean;
      vv = sopt.PCinfo.direction;
      dd = sopt.PCinfo.variance;
   end
   dcum = cumsum(dd);
   nPC = find(dcum >= p1,1);
else
   nPC = nVar;
   vv = eye(nVar);
end
vPC = vv(:,1:nPC);
if ~isempty(Dinfo.D)
   if size(Dinfo.D,1) == nPC
      A = Dinfo.D;
      uq = Dinfo.uq;
   else
      A = [];
      uq = [];
   end
end
if nPC < nVar
   if strcmp(sopt.UncertaintyEstimation,'Outer')
      sopt.UncertaintyEstimation = 'Inner';
   end
   uqv = zeros(nPC,2);
   switch sopt.UncertaintyEstimation
      case 'Inner'
         opt.Prediction = 'inner';
         for i = 1:nPC
            tv = [-xAve*vv(:,i); vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         for i = 1:nPC
            f = xC * vv(:,i);
            uqv(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         sf = 1-sopt.UncertaintyTruncation;
         for i = 1:nPC
            f = xC * vv(:,i);
            uqv(i,:) = sf*[min(f) max(f)];
         end
   end
   sVar = generateVar([],uqv);
else
   sVar = obj.Variables;
end
nh = round(0.04*nCut);
uq2 = zeros(nCut,2);
D = randn(nPC,nCut);
D = normc(D);
switch sopt.UncertaintyEstimation
   case 'Outer'
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nCut
         tv = [0; D(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq2(i,:) = [preQ.min preQ.max];
         if mod(i,nh) == 1
            waitbar(i/nCut,h);
         end
      end
      close(h);
   case 'Inner'
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nCut
         if nPC < nVar
            uq2(i,:) = obj.predictDirection_sample(vv,nPC,xAve,D(:,i));
         else
            tv = [0; D(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
         if mod(i,nh) == 1
            waitbar(i/nCut,h);
         end
      end
      close(h);
   case 'Sample'
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nCut
         tv = vPC*D(:,i);
         if nPC < nVar
            f = xC*tv;
         else
            f = xS*tv;
         end
         uq2(i,:) = [min(f) max(f)];
         if mod(i,nh) == 1
            waitbar(i/nCut,h);
         end
      end
      close(h);
   case 'Truncation'
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nCut
         tv = vPC*D(:,i);
         f = xC*tv;
         uq2(i,:) = sf*[min(f) max(f)];
         if nPC == nVar
            uq2(i,:) = uq2(i,:) + xAve*tv;
         end
         if mod(i,nh) == 1
            waitbar(i/nCut,h);
         end
      end
      close(h);
end
% t1 = toc;
A = [A D]';
uq = [uq;uq2];
if size(A,1) > 0
   sVar = sVar.addLinearConstraint(A,uq(:,1),uq(:,2));
end
if sopt.ParameterScaling
   sVar = sVar.calScale;
end
% t2 =toc;
% tt = t2-t1;
ic = 1;
n1 = size(xVal,1);
ns = sopt.BatchMaxSample;
h = waitbar(n1/N,'Collecting uniform samples of the feasible set...');
while n1 < N
   xcand = sVar.collectSamples(ns,[],opt.SampleOption);
   if nPC < nVar
      xtmp = repmat(xAve*vv,ns,1);
      xtmp(:,1:nPC) = xtmp(:,1:nPC) + xcand;
      xtmp = xtmp*vv';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xVal = [xVal; xtmp(iF,:)];
   n1 = size(xVal,1);
   waitbar(min(n1/N,1),h);
   if n1 >= N
      eff = n1/ic/ns;
      xVal = xVal(randperm(n1,N),:);
      close(h);
      return;
   end
   if ic == 5
      if n1 <= th*ic*ns
         eff = n1/ic/ns;
         disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
         close(h)
         newCut = input('How many extra cut do you want to include? \n');
         if newCut >= 0
            Dinfo.D = A';
            Dinfo.uq = uq;
            opt.SampleOption.ExtraCut = ceil(newCut);
            [xVal2,eff,opt] = RM_hr(obj,N-n1,xS,Dinfo,opt);
            xVal = [xVal; xVal2];
            xVal = xVal(randperm(size(xVal,1),size(xVal,1)),:);
            return;
         else
            return;
         end
      end
   end
end

