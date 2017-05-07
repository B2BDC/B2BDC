function [xVal,eff,opt] = RM_sobol_par(obj,N,xF,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 3, 2016     Wenyu Li

xVal = [];
sOpt = opt.SampleOption;
th2 = sOpt.RejectionTol;
s1 = sOpt.UncertaintyEstimation;
p1 = 1-sOpt.PCTruncation;
vars = obj.Variables;
nVar = vars.Length;
% t1 = toc;
if isempty(xF) && (isempty(sOpt.PCinfo) || strcmp(s1,'Sample') || strcmp(s1,'Truncation'))
   nS = min(10^7,10*nVar^3);
   xF = obj.collectSamples(nS,[],opt);
else
   nS = size(xF,1);
end
if ~isempty(xF)
   xAve = mean(xF);
   xC = xF - repmat(xAve,nS,1);
   [V,D] = eig((xC' * xC)/nS);
   Ddiag = diag(D);
   [dd,id] = sort(Ddiag,'descend');
   dd = dd/sum(dd);
   vv = V(:,id);
   pcinfo.direction = vv;
   pcinfo.variance = dd;
   pcinfo.mean = xAve;
   opt.SampleOption.PCinfo = pcinfo;
else
   xAve = sOpt.PCinfo.mean;
   vv = sOpt.PCinfo.direction;
   dd = sOpt.PCinfo.variance;
end
dcum = cumsum(dd);
nPC = find(dcum >= p1,1);
uq = zeros(nPC,2);
p = parpool('IdleTimeout', 120);
npool = p.NumWorkers;
if npool >= 4
   nt1 = npool;
else
   nt1 = 2*npool;
end
disp('Calculating directional uncertainty...');
switch s1
   case 'Outer'
      opt.Prediction = 'outer';
      opt.ExtraLinFraction = -1;
      parfor i = 1:nPC
         tv = [-xAve*vv(:,i); vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq(i,:) = [preQ.min, preQ.max];
      end
      close(h);
   case 'Inner'
      opt.Prediction = 'inner';
      parfor i = 1:nPC
         tv = [-xAve*vv(:,i); vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq(i,:) = [preQ.min, preQ.max];
      end
      close(h);
   case 'Sample'
      parfor i = 1:nPC
         f = xC*vv(:,i);
         uq(i,:) = [min(f), max(f)];
      end
      close(h);
   case 'Truncation'
      sf = 1-sOpt.UncertaintyTruncation;
      parfor i = 1:nPC
         f = xC*vv(:,i);
         uq(i,1) = [sf*min(f), sf*max(f)];
      end
end
% t2 = toc;
% tt = t2-t1;
ns = sOpt.BatchMaxSample;
xmin = uq(:,1)';
dx = uq(:,2)-uq(:,1);
tic;
xCand = cell(nt1,1);
disp('Gathering samples...');
parfor i = 1:nt1
   sob = sobolset(nPC,'Skip',round(ns*100*rand(1)),'Leap',100);
   sob = scramble(sob,'MatousekAffineOwen');
   rr = net(sob,ns);
   xtmp = zeros(ns,nVar);
   xtmp(:,1:nPC) = repmat(xmin,ns,1)+rr.*repmat(dx',ns,1);
   xtmp = xtmp*vv'+repmat(xAve,ns,1);
   iF = obj.isFeasiblePoint(xtmp);
   xCand{i} = xtmp(iF,:);
end
for i = 1:nt1
   xVal = [xVal; xCand{i}];
end
% t3 = toc;
n1 = size(xVal,1);
if n1 <= th2*nt1*ns
   eff = n1/nt1/ns;
   disp(['Numerical efficiency of current sampling method is' num2str(eff)])
   delete(p)
   return    
else
   if n1 >= N
      eff = n1/nt1/ns;
      xVal = xVal(randperm(n1,N),:);
      return
   end
   eff = n1/ns/nt1;
   nEst = ceil(1.2*(N-n1)/eff/ns);
   disp(['The estimated batch of runs is: ' num2str(nEst)])
   disp(['The estimated running time is: ' num2str(t3*nEst/nt1/60) 'minutes'])
end
xCand = cell(nEst,1);
parfor i = 1:nEst
   sob = sobolset(nPC,'Skip',round(ns*10*nEst*rand(1)),'Leap',100);
   sob = scramble(sob,'MatousekAffineOwen');
   rr = net(sob,ns);
   xtmp = zeros(ns,nVar);
   xtmp(:,1:nPC) = repmat(xmin,ns,1)+rr.*repmat(dx',ns,1);
   xtmp = xtmp*vv'+repmat(xAve,ns,1);
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
      sob = sobolset(nPC,'Skip',round(ns*20*nEst*rand(1)),'Leap',100);
      sob = scramble(sob,'MatousekAffineOwen');
      rr = net(sob,ns);
      xtmp = zeros(ns,nVar);
      xtmp(:,1:nPC) = repmat(xmin,ns,1)+rr.*repmat(dx',ns,1);
      xtmp = xtmp*vv'+repmat(xAve,ns,1);
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
   for i = 1:length(xCand)
      xVal = [xVal; xCand{i}];
   end
else
   is = 0;
end
delete(p)
n1 = size(xVal,1);
itt = (nt1+nEst+is*ceil(0.2*nEst))*ns;
if n1 < N
   disp('Not enough samples are found within 120% calculating budget')
end
eff = n1/itt;
xVal = xVal(randperm(n1,min(n1,N)),:);