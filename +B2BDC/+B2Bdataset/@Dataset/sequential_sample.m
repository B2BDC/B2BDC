function [xVal,eff] = sequential_sample(obj,nPC,xF,Dinfo,opt)

% [XVAL,EFF] = SEQUENTIAL_SAMPLE(OBJ,N,XS,) generates uniform samples by
% the sequential sampling strategy.

%  Created: April 23, 2018     Wenyu Li

if ~obj.isConsistent(opt)
   error('The dataset is inconsistent');
end
sopt = opt.SampleOption;
vars = obj.Variables;
nVar = vars.Length;
s1 = sopt.UncertaintyEstimation;
da = sopt.DirectionAdaption;
nCut = sopt.ExtraCut;
p1 = sopt.DataStorePath;
if isempty(p1)
   disp('The path of previous samples is not specified');
   xVal = [];
   eff = 0;
   return
else
   p2 = fullfile(p1,['Sequential sample, nPC = ' num2str(nPC(1)) '+' num2str(nPC(2))]);
end
if da > 0
   daFlag = true;
else
   daFlag = false;
end
fid = fopen(fullfile(p2,'Computation performance.txt'),'w');
if (strcmp(s1,'Sample') || strcmp(s1,'Truncation')) && isempty(xF)
   nF = min(10^6,10*nVar^3);
   xF = obj.collectSamples(nF,[],opt);
else
   nF = size(xF,1);
end
rr = load(fullfile(p1,'xSample'));
xS = rr.xVal;
nS = size(xS,1);
rr = load(fullfile(p1,'PCinfo'));
PCinfo = rr.sampledPC;
V = PCinfo.V;
xAve = PCinfo.x0;
xC = xF-repmat(xAve,nF,1);
vv = V(:,nPC(1)+1:sum(nPC));
D = randn(nPC(2),nCut);
D = normc(D);
uqs = zeros(nCut,2);
uqv = zeros(nPC(2),:);
opt.Display = false;
tic;
switch s1
   case 'Outer'
      opt.Prediction = 'outer';
      for i = 1:nPC(2)
         tv = [-xAve*vv(:,i);vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uqv(i,:) = [preQ.min preQ.max];
      end
      for i = 1:nCut
         tv = [-xAve*vv*D(:,i); vv*D(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uqs(i,:) = [preQ.min preQ.max];
      end
   case 'Inner'
      opt.Prediction = 'inner';
      for i = 1:nPC(2)
         tv = [-xAve*vv(:,i);vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uqv(i,:) = [preQ.min preQ.max];
      end
      for i = 1:nCut
         tv = [-xAve*vv*D(:,i); vv*D(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uqs(i,:) = [preQ.min preQ.max];
      end
   case 'Sample'
      for i = 1:nPC(2)
         f = xC*vv(:,i);
         uqv(i,:) = [min(f) max(f)];
      end
      for i = 1:nCut
         f = xC*vv*D(:,i);
         uqs(i,:) = [min(f) max(f)];
      end
   case 'Truncation'
      opt.Prediction = 'inner';
      for i = 1:nPC(2)
         tv = [-xAve*vv(:,i);vv(:,i)];
         tM = generateModel(tv,vars);
         preQ = obj.predictQOI(tM,opt);
         uqv(i,:) = [preQ.min preQ.max];
      end
      for i = 1:nCut
         tv = [-xAve*vv*D(:,i); vv*D(:,i)];
         tM = generateModel(tv,vars);
         preQ = subDS.predictQOI(tM,opt);
         fmin = preQ.min;
         fmax = preQ.max;
         uqs(i,:) = sf*[fmin, fmax];
      end
end
tt = toc;
fprintf(fid,'The CPU time to calcualte the polytope with %d directions is: %5.3E \n',nCut,tt);
sVar = generateVar([],uqv);
sVar = sVar.addLinearConstraint(D',uqs(:,1),uqs(:,2));
nT = round(0.05*nS);
tic;
while any(daFlag)
   nBatch = sopt.BatchMaxSample;
   idT = randperm(nS,nT);
   xcand = [];
   neff = 0;
   for i = idT
      if isempty(xcand)
         xcand = sVar.collectSamples(nBatch,[],opt.SampleOption);
      else
         xcand = sVar.collectSamples(nBatch,xcand(end,:),opt.SampleOption);
      end
      xtmp = repmat(xS(i,:),nBatch,1) + xcand*vv';
      flag = obj.isFeasiblePoint(xtmp);
      neff = neff+sum(flag);
   end
   ss = neff/nT;
   if ss<0.7
      error('The sampling efficiency is too low')
   elseif ss>1.2
      nBatch = 1.2*nBatch/ss;
   end
   xApp = [];
   nda = round(da*nS);
   idT = randperm(nS,nda);
   for i = idT
      if isempty(xcand)
         xcand = sVar.collectSamples(nBatch,[],opt.SampleOption);
      else
         xcand = sVar.collectSamples(nBatch,xcand(end,:),opt.SampleOption);
      end
      xtmp = repmat(xS(i,:),nBatch,1) + xcand*vv';
      flag = obj.isFeasiblePoint(xtmp);
      xApp = [xApp; xtmp(flag,:)];
   end
   [sVar,daFlag] = sVar.DirectionalAdapt(xApp,daFlag);
end
nBatch = sopt.BatchMaxSample;
idT = randperm(nS,nT);
xcand = [];
neff = 0;
for i = idT
   if isempty(xcand)
      xcand = sVar.collectSamples(nBatch,[],opt.SampleOption);
   else
      xcand = sVar.collectSamples(nBatch,xcand(end,:),opt.SampleOption);
   end
   xtmp = repmat(xS(i,:),nBatch,1) + xcand*vv';
   flag = obj.isFeasiblePoint(xtmp);
   neff = neff+sum(flag);
end
ss = neff/nT;
if ss<0.7
   error('The sampling efficiency is too low')
elseif ss>1.2
   nBatch = 1.2*nBatch/ss;
end
tt = toc;
fprintf(fid,'The CPU time to adapt sample batch size and directional facets is: %5.3E \n',tt);
xVal = [];
nTime = round(0.05*nS);
tic;
for i = 1:nS
   if isempty(xcand)
      xcand = sVar.collectSamples(nBatch,[],opt.SampleOption);
   else
      xcand = sVar.collectSamples(nBatch,xcand(end,:),opt.SampleOption);
   end
   xtmp = repmat(xS(i,:),nBatch,1) + xcand*vv';
   flag = obj.isFeasiblePoint(xtmp);
   xVal = [xVal; xtmp(flag,:)];
   if mod(i,nTime) == 0
      tt = toc;
      fprintf('Current progress: %4.2f%%, estimated time left: %5.3E \n',i/nS,tt*nS/i);
   end
end
save(fullfile(p2,'xSample'),'xVal');
eff = size(xVal,1)/nS;
save(fullfile(p2,'efficiency'),'eff');