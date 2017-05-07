function [xVal,eff,opt] = RM_hr_par(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 18, 2016     Wenyu Li

% tic;
th = opt.SampleOption.RejectionTol;
p1 = opt.SampleOption.PCTruncation;
xVal = [];
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
sVar = vars;
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
      vv = V(:,id);
      pcinfo.direction = vv;
      pcinfo.variance = dd;
      pcinfo.mean = xAve;
      opt.SampleOption.PCinfo = pcinfo;
   else
      xAve = sopt.PCinfo.mean;
      vv = sopt.PCinfo.direction;
      dd = sopt.PCinfo.variance;
   end
   dcum = cumsum(dd);
   dcum = dcum/dcum(end);
   nPC = find(dcum >= 1 - p1,1);
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
else
   A = [];
   uq = [];
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
uq2 = zeros(nCut,2);
D = randn(nPC,nCut);
D = normc(D);
p = parpool('IdleTimeout', 120);
npool = p.NumWorkers;
if npool >= 4
   nt1 = npool;
else
   nt1 = 2*npool;
end
% t1 = toc;
switch sopt.UncertaintyEstimation
   case 'Outer'
      parfor i = 1:nCut
         tv = [0; D(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq2(i,:) = [preQ.min preQ.max];
      end
   case 'Inner'
      parfor i = 1:nCut
         if nPC < nVar
            uq2(i,:) = obj.predictDirection_sample(vv,nPC,xAve,D(:,i));
         else
            tv = [0; D(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      end
   case 'Sample'
      parfor i = 1:nCut
         tv = vPC*D(:,i);
         if nPC < nVar
            f = xC*tv;
         else
            f = xS*tv;
         end
         uq2(i,:) = [min(f) max(f)];
      end
   case 'Truncation'
      parfor i = 1:nCut
         tv = vPC*D(:,i);
         f = xC*tv;
         uq2(i,:) = sf*[min(f) max(f)];
         if nPC == nVar
            uq2(i,:) = uq2(i,:) + xAve*tv;
         end
      end
end
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
ns = sopt.BatchMaxSample;
tic;
xCand = cell(nt1,1);
parfor i = 1:nt1
   xcand = sVar.collectSamples_par(ns,[],opt.SampleOption);
   if nPC < nVar
      xtmp = repmat(xAve*vv,ns,1);
      xtmp(:,1:nPC) = xtmp(:,1:nPC) + xcand;
      xtmp = xtmp*vv';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xCand{i} = xtmp(iF,:);
end
for i = 1:nt1
   xVal = [xVal; xCand{i}];
end
% t3 = toc;
n1 = size(xVal,1);
if n1 <= th*nt1*ns
   eff = n1/nt1/ns;
   disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
   delete(p)
   return;
else
   if n1 >= N
      xVal = xVal(randperm(n1,N),:);
      eff = n1/nt1/ns;
      delete(p)
      return;
   end
   eff = n1/ns/nt1;
   nEst = ceil(1.2*(N-n1)/eff/ns);
   disp(['The estimated batch of runs is: ' num2str(nEst)])
   disp(['The estimated running time is: ' num2str(t3*nEst/60/nt1) 'minutes'])
end
xCand = cell(nEst,1);
parfor i = 1:nEst
   xcand = sVar.collectSamples_par(ns,[],opt.SampleOption);
   if nPC < nVar
      xtmp = repmat(xAve*vv,ns,1);
      xtmp(:,1:nPC) = xtmp(:,1:nPC) + xcand;
      xtmp = xtmp*vv';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xCand{i} = xtmp(iF,:);
end
for i = 1:nEst
   xVal = [xVal; xCand{i}];
end
n1 = size(xVal,1);
if n1 < N
   is = 1;
   xCand = cell(ceil(0.2*nEst),1);
   parfor i = 1:length(xCand)
      xcand = sVar.collectSamples_par(ns,[],opt.SampleOption);
      if nPC < nVar
         xtmp = repmat(xAve*vv,ns,1);
         xtmp(:,1:nPC) = xtmp(:,1:nPC) + xcand;
         xtmp = xtmp*vv';
      else
         xtmp = xcand;
      end
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
   for i = 1:length(xCand)
      xVal = [xVal; xCand{i}];
   end
else
   is = 0;
end
delete(p);
n1 = size(xVal,1);
itt = (nt1+nEst+is*ceil(0.2*nEst))*ns;
if n1 < N
   disp('Not enough samples are found within 120% calculating budget')
end
eff = n1/itt;
xVal = xVal(randperm(n1,min(n1,N)),:);
