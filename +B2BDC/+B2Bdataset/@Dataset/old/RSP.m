function [xVal,eff] = RSP(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 3, 2016     Wenyu Li
%  Modified Nov 22, 2017  (save result during the sample)  Wenyu Li
%  Modified April 23, 2018 (adaptive update of polytope)  Wenyu Li

th = opt.SampleOption.RejectionTol;
da = opt.SampleOption.DirectionAdaption;
if da > 0
   daFlag = true;
else
   daFlag = false;
end
xVal = [];
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nPC = nVar-opt.SampleOption.TruncatedPC;
if ~isempty(sopt.DataStorePath) && isdir(sopt.DataStorePath)
   savePath = sopt.DataStorePath;
   fid = fopen(fullfile(savePath,'Computation performance.txt'),'w');
else
   savePath = [];
   fid = -1;
end
if nPC < 1
   nPC = 1;
end
nCut = sopt.ExtraCut;
if nPC < nVar  ||...
      ((strcmp(s1,'Sample') || strcmp(s1,'Truncation')) && nCut >0)
   if isempty(xS)
      nS = min(10^6,10*nVar^3);
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
   [V,DD] = eig(cov(xS));
   Ddiag = diag(DD);
   [~,id] = sort(Ddiag,'descend');
   vv = V(:,id);
%    y0 = xAve*vv;
   sampledPC.V = vv;
   sampledPC.D = Ddiag(id);
   sampledPC.x0 = xAve;
   if ~isempty(savePath)
      save(fullfile(savePath,'PCinfo'),'sampledPC');
   end
else
   vv = eye(nVar);
end
if ~isempty(Dinfo)
   if size(Dinfo.D,1) == nPC
      A = Dinfo.D;
      uq = Dinfo.uq;
   else
      disp('The input D matrix should have nPC rows');
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
   projectDS = obj.projectDSonActiveSubspace(sampledPC,nPC);
   vars = projectDS.Variables;
   uqv = zeros(nPC,2);
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         for i = 1:nPC
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         for i = 1:nPC
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
%          xC = projectDS.collectSamples(1e5+1e5,[],opt);
%          xC = xC(1e5+1:end,:);
         xx = xS-repmat(mean(xS),size(xS,1),1);
         xC = xx*vv(:,1:nPC);
%          if ~isempty(savePath)
%             save(fullfile(savePath,'xF_pc'),'xC');
%          end
         uqv = [min(xC)' max(xC)'];
      case 'Truncation'
         sf = 1-sopt.UncertaintyTruncation;
         opt.Prediction = 'inner';
         for i = 1:nCut
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            fmin = preQ.min;
            fmax = preQ.max;
            uqv(i,:) = sf*[fmin, fmax];
         end
%          for i = 1:nPC
%             f = xC * vv(:,i);
%             f0 = xAve * vv(:,i);
%             uqv(i,:) = sf*[min(f) max(f)];
%          end
   end
   sVar = generateVar([],uqv);
else
   sVar = obj.Variables;
   if strcmp(sopt.UncertaintyEstimation,'Sample')
      xC = xS;
   end
end
if nCut > 0
   tic;
   if nPC < nVar
      subDS = projectDS;
   else
      subDS = obj;
   end
   uq2 = zeros(nCut,2);
   D = randn(nPC,nCut);
   D = normc(D);
   opt.Display = false;
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         for i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         for i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         for i = 1:nCut
            f = xC*D(:,i);
            uq2(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         opt.Prediction = 'inner';
         for i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            fmin = preQ.min;
            fmax = preQ.max;
            uq2(i,:) = sf*[fmin, fmax];
         end
%          xC_pc = xC*vv(:,1:nPC);
%          xAve_pc = mean(xS_pc);
%          for i = 1:nCut
%             f = xC_pc*D(:,i);
%             mf = mean(f);
%             mf = xAve_pc*D(:,i);
%             uq2(i,:) = [mf+(1-sf)*(min(f)-mf), mf+(1-sf)*(max(f)-mf)];
%          end
   end
   tt = toc;
   if fid > 0
      fprintf(fid,'The CPU time to calcualte the polytope with %d directions is: %5.3E \n',nCut,tt);
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
      save(fullfile(savePath,'polytopeInfo'),'sVar');
   end
end
% if sopt.ParameterScaling
%    sVar = sVar.calScale;
% end
ns = sopt.BatchMaxSample;
tic;
while any(daFlag)
   ic = 0;
   xApp = [];
   xcand = [];
   n1 = size(xApp,1);
   while n1 < da*N
      if isempty(xcand)
         xcand = sVar.collectSamples(ns,[],opt.SampleOption);
      else
         xcand = sVar.collectSamples(ns,xcand(end,:),opt.SampleOption);
      end
      ic = ic+1;
      if nPC < nVar
         xtmp = repmat(xAve,ns,1);
         xtmp = xtmp+xcand*vv(:,1:nPC)';
      else
         xtmp = xcand;
      end
      iF = obj.isFeasiblePoint(xtmp);
      xApp = [xApp; xtmp(iF,:)];
      if ic == 5
         if n1 <= th*ic*ns
            eff = n1/ic/ns;
            disp(['Numerical efficiency of current sampling method is ' num2str(eff)]);
            return;
         end
      end
   end
   [sVar,daFlag] = sVar.DirectionalAdapt(xApp,daFlag);
end
tt = toc;
if fid > 0
   fprintf(fid,'The CPU time to update polytope directional facets is: %5.3E \n',tt);
end
tic;
ic = 0;
xcand = [];
xVal = [];
n1 = size(xVal,1);
while n1 < N
   if isempty(xcand)
      xcand = sVar.collectSamples(ns,[],opt.SampleOption);
   else
      xcand = sVar.collectSamples(ns,xcand(end,:),opt.SampleOption);
   end
   ic = ic+1;
   if nPC < nVar
      xtmp = repmat(xAve,ns,1);
      xtmp = xtmp+xcand*vv(:,1:nPC)';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xVal = [xVal; xtmp(iF,:)];
   if ~isempty(savePath)
      save(fullfile(savePath,'xSample'),'xVal');
   end
   n1 = size(xVal,1);
   if n1 >= N
      eff = n1/ic/ns;
      if ~isempty(savePath)
         save(fullfile(savePath,'efficiency'),'eff');
      end
      xVal = xVal(randperm(n1,N),:);
      xVal = xVal(1:N,:);
      if fid > 0
         fclose(fid);
      end
      return;
   end
   if ic == 5
      if n1 <= th*ic*ns
         eff = n1/ic/ns;
         disp(['Numerical efficiency of current sampling method is ' num2str(eff)]);
%          newCut = input('How many extra cut do you want to include? \n (negative value to exit) \n' );
%          if newCut > 0
%             Dinfo.D = A';
%             Dinfo.uq = uq;
%             opt.SampleOption.ExtraCut = ceil(newCut);
%             [xVal2,eff] = RSP(obj,N-n1,xS,Dinfo,opt);
%             xVal = [xVal; xVal2];
%             xVal = xVal(randperm(size(xVal,1),size(xVal,1)),:);
         if fid > 0
            fclose(fid);
         end
         return;
      end
   end
   tt = toc;
   if fid > 0
      fprintf('Current progress: %4.2f%%, estimated time left: %5.3E \n',n1/N,tt*N/n1);
   end
end

