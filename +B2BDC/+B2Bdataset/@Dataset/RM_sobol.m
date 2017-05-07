function [xVal,eff,opt] = RM_sobol(obj,N,xF,opt)
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
% nPC = sum(dd >= p1*dd(1));
uq = zeros(nPC,2);
switch s1
   case 'Outer'
      opt.Prediction = 'outer';
      opt.ExtraLinFraction = -1;
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nPC
         tv = [-xAve*vv(:,i); vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq(i,1) = preQ.min;
         uq(i,2) = preQ.max;
         waitbar(i/nPC,h);
      end
      close(h);
   case 'Inner'
      opt.Prediction = 'inner';
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nPC
         tv = [-xAve*vv(:,i); vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uq(i,1) = preQ.min;
         uq(i,2) = preQ.max;
         waitbar(i/nPC,h);
      end
      close(h);
   case 'Sample'
      h = waitbar(0,'Calculating directional uncertainty...');
      for i = 1:nPC
         f = xC*vv(:,i);
         uq(i,1) = min(f);
         uq(i,2) = max(f);
         waitbar(i/nPC,h);
      end
      close(h);
   case 'Truncation'
      h = waitbar(0,'Calculating directional uncertainty...');
      sf = 1-sOpt.UncertaintyTruncation;
      for i = 1:nPC
         f = xC*vv(:,i);
         uq(i,1) = sf*min(f);
         uq(i,2) = sf*max(f);
         waitbar(i/nPC,h);
      end
      close(h);
end
% t2 = toc;
% tt = t2-t1;
ic = 1;
n1 = size(xVal,1);
ns = sOpt.BatchMaxSample;
h = waitbar(n1/N,'Collecting uniform samples of the feasible set...');
while n1 < N
   xmin = uq(:,1)';
   dx = uq(:,2)-uq(:,1);
   sob = sobolset(nPC,'Skip',ns*ic,'Leap',100);
   sob = scramble(sob,'MatousekAffineOwen');
   rr = net(sob,ns);
   xtmp = zeros(ns,nVar);
   xtmp(:,1:nPC) = repmat(xmin,ns,1)+rr.*repmat(dx',ns,1);
   xtmp = xtmp*vv'+repmat(xAve,ns,1);
   iF = obj.isFeasiblePoint(xtmp);
   xVal = [xVal; xtmp(iF,:)];
   n1 = size(xVal,1);
   waitbar(min(n1/N,1),h);
   if n1 >= N
      eff = n1/ic/ns;
      close(h);
      break
   end
   if ic == 5
      if n1 <= th2*ic*ns
         eff = n1/ic/ns;
         disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
         return
      end
   end
   ic = ic+1;
end
xVal = xVal(randperm(n1,N),:);