function [xVal,eff] = RSP_par(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 18, 2016     Wenyu Li


p = parpool('IdleTimeout', 120);
npool = p.NumWorkers;
if npool >= 4
   nt1 = npool;
else
   nt1 = 2*npool;
end
th = opt.SampleOption.RejectionTol;
xVal = [];
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nCut = sopt.ExtraCut;
nPC = nVar-opt.SampleOption.TruncatedPC;
if ~isempty(sopt.DataStorePath) && isdir(sopt.DataStorePath)
   savePath = sopt.DataStorePath;
else
   savePath = [];
end
if nPC < 1
   nPC = 1;
end
if nPC < nVar ||...
      ((strcmp(s1,'Sample') || strcmp(s1,'Truncation')) && nCut >0)
   if isempty(xS)
      nS = min(10^7,10*nVar^3);
      xS = obj.collectSamples(nS,[],opt);
   else
      nS = size(xS,1);
   end
end
if nPC < nVar
   xAve = mean(xS);
   if ~obj.isFeasiblePoint(xAve)
      Xdd = xS - repmat(xAve,nS,1);
      xdd = sqrt(sum(Xdd.^2,2));
      [~,imin] = min(xdd);
      xAve = xS(imin,:);
   end
   xC = xS - repmat(xAve,nS,1);
   [V,DD] = eig((xC' * xC)/nS);
   Ddiag = diag(DD);
   [~,id] = sort(Ddiag,'descend');
   vv = V(:,id);
   y0 = xAve*vv;
   if ~isempty(savePath)
      sampledPC.V = V;
      sampledPC.D = diag(Ddiag);
      save(fullfile(savePath,'PCinfo'),'sampledPC');
   end
else
   vv = eye(nVar);
end
if ~isempty(Dinfo)
   if size(Dinfo.D,1) == nPC
      disp('The input D matrix should have nPC rows');
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
%    if strcmp(sopt.UncertaintyEstimation,'Outer')
%       sopt.UncertaintyEstimation = 'Inner';
%    end
   uqv = zeros(nPC,2);
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         parfor i = 1:nPC
            tv = [-y0(i); vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         parfor i = 1:nPC
            tv = [-y0(i); vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         parfor i = 1:nPC
            f = xC * vv(:,i);
            uqv(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         sf = 1-sopt.UncertaintyTruncation;
         parfor i = 1:nPC
            f = xC * vv(:,i);
%             f0 = xAve * vv(:,i);
            uqv(i,:) = sf*[min(f) max(f)];
         end
   end
   sVar = generateVar([],uqv);
else
   sVar = obj.Variables;
end
if nCut > 0
   if nPC < nVar
      subDS = generateSubDS(obj,vv,nPC,xAve,sVar);
   else
      subDS = obj;
   end
   uq2 = zeros(nCut,2);
   D = randn(nPC,nCut);
   D = normc(D);
   % t1 = toc;
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         parfor i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         parfor i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         xC_pc = xC*vv(:,1:nPC);
         parfor i = 1:nCut
            f = xC_pc*D(:,i);
            uq2(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         xC_pc = xC*vv(:,1:nPC);
%          xAve_pc = mean(xS_pc);
         parfor i = 1:nCut
            f = xC_pc*D(:,i);
            mf = mean(f);
            uq2(i,:) = [mf+(1-sf)*(min(f)-mf) mf+(1-sf)*(max(f)-mf)];
         end
   end
else
   D = [];
   uq2 = [];
end
A = [A D]';
uq = [uq;uq2];
if size(A,1) > 0
   sVar = sVar.addLinearConstraint(A,uq(:,1),uq(:,2));
   if ~isempty(savePath)
      polytope.direction = A;
      polytope.bounds = uq;
      save(fullfile(savePath,'polytopeInfo'),'polytope');
   end
end
% if sopt.ParameterScaling
%    sVar = sVar.calScale;
% end
ns = sopt.BatchMaxSample;
xCand = cell(nt1,1);
if nPC < nVar
   parfor i = 1:nt1
      xcand = sVar.collectSamples(ns,[],opt.SampleOption);
      xtmp = repmat(y0,ns,1);
      xtmp(:,1:nPC) = xtmp(:,1:nPC)+xcand;
      xtmp = xtmp*vv';
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
else
   parfor i = 1:nt1
      xtmp = sVar.collectSamples(ns,[],opt.SampleOption);
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
end
for i = 1:nt1
   xVal = [xVal; xCand{i}];
end
if ~isempty(savePath)
   save(fullfile(savePath,'xSample'),'xVal');
end
n1 = size(xVal,1);
if n1 <= th*nt1*ns
   eff = n1/nt1/ns;
   disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
   delete(p)
   return;
else
   if n1 >= N
      xVal = xVal(1:N,:);
      eff = n1/nt1/ns;
      if ~isempty(savePath)
         save(fullfile(savePath,'efficiency'),'eff');
      end
      delete(p)
      return;
   end
   eff = n1/ns/nt1;
   if ~isempty(savePath)
      save(fullfile(savePath,'efficiency'),'eff');
   end
   nEst = ceil(1.2*(N-n1)/eff/ns);
   disp(['The estimated batch of runs is: ' num2str(nEst)])
   disp(['The estimated running time is: ' num2str(t3*nEst/60/nt1) 'minutes'])
end
xCand = cell(nEst,1);
if nPC < nVar
   parfor i = 1:nEst
      xcand = sVar.collectSamples(ns,[],opt.SampleOption);
      xtmp = repmat(y0,ns,1);
      xtmp(:,1:nPC) = xtmp(:,1:nPC)+xcand;
      xtmp = xtmp*vv';
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
else
   parfor i = 1:nEst
   xtmp = sVar.collectSamples(ns,[],opt.SampleOption);
   iF = obj.isFeasiblePoint(xtmp);
   xCand{i} = xtmp(iF,:);
   end
end
for i = 1:nEst
   xVal = [xVal; xCand{i}];
end
if ~isempty(savePath)
   save(fullfile(savePath,'xSample'),'xVal');
end
n1 = size(xVal,1);
if n1 < N
   is = 1;
   xCand = cell(ceil(0.2*nEst),1);
   if nPC < nVar
      parfor i = 1:length(xCand)
         xcand = sVar.collectSamples(ns,[],opt.SampleOption);
         xtmp = repmat(y0,ns,1);
         xtmp(:,1:nPC) = xtmp(:,1:nPC)+xcand;
         xtmp = xtmp*vv';
         iF = obj.isFeasiblePoint(xtmp);
         xCand{i} = xtmp(iF,:);
      end
   else
      parfor i = 1:length(xCand)
         xtmp = sVar.collectSamples(ns,[],opt.SampleOption);
         iF = obj.isFeasiblePoint(xtmp);
         xCand{i} = xtmp(iF,:);
      end
   end
   for i = 1:length(xCand)
      xVal = [xVal; xCand{i}];
   end
   if ~isempty(savePath)
      save(fullfile(savePath,'xSample'),'xVal');
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
if ~isempty(savePath)
   save(fullfile(savePath,'efficiency'),'eff');
end
xVal = xVal(1:min(n1,N),:);
