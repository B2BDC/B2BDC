function [xVal,eff] = RSP(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 3, 2016     Wenyu Li
%  Modified Nov 22, 2017  (save result during the sample)  Wenyu Li

% tic;
th = opt.SampleOption.RejectionTol;
xVal = [];
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nPC = nVar-opt.SampleOption.TruncatedPC;
if ~isempty(sopt.DataStorePath) && isdir(sopt.DataStorePath)
   savePath = sopt.DataStorePath;
else
   savePath = [];
end
if nPC < 1
   nPC = 1;
end
nCut = sopt.ExtraCut;
if nPC < nVar  ||...
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
   xC = xS - repmat(xAve,nS,1);
   [V,DD] = eig((xC' * xC)/nS);
   Ddiag = diag(DD);
   [~,id] = sort(Ddiag,'descend');
   vv = V(:,id);
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
         for i = 1:nPC
            tv = [0; vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         for i = 1:nPC
            tv = [0; vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         for i = 1:nPC
            f = xS * vv(:,i);
            uqv(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         sf = 1-sopt.UncertaintyTruncation;
         for i = 1:nPC
            f = xC * vv(:,i);
            f0 = xAve * vv(:,i);
            uqv(i,:) = sf*[min(f) max(f)] + f0;
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
         xS_pc = xS*vv(:,1:nPC);
         for i = 1:nCut
            f = xS_pc*D(:,i);
            uq2(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         xS_pc = xS*vv(:,1:nPC);
         xAve_pc = mean(xS_pc);
         for i = 1:nCut
            f = xS*D(:,i);
            mf = xAve_pc*D(:,i);
            uq2(i,:) = [mf+(1-sf)*(min(f)-mf) mf+(1-sf)*(max(f)-mf)];
         end
   end
else
   D = [];
   uq2 = [];
end
% t1 = toc;
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
% t2 =toc;
% tt = t2-t1;
ic = 1;
n1 = size(xVal,1);
ns = sopt.BatchMaxSample;
h = waitbar(n1/N,'Collecting uniform samples of the feasible set...');
while n1 < N
   xcand = sVar.collectSamples(ns,[],opt.SampleOption);
   if nPC < nVar
      xtmp = [zeros(ns,nPC) repmat(xAve*vv(:,nPC+1:end),ns,1)];
      xtmp(:,1:nPC) = xcand;
      xtmp = xtmp*vv';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xVal = [xVal; xtmp(iF,:)];
   if ~isempty(savePath)
      save(fullfile(savePath,'xSample'),'xVal');
   end
   n1 = size(xVal,1);
   waitbar(min(n1/N,1),h);
   if n1 >= N
      eff = n1/ic/ns;
      if ~isempty(savePath)
         save(fullfile(savePath,'efficiency'),'eff');
      end
      xVal = xVal(randperm(n1,N),:);
      xVal = xVal(1:N,:);
      close(h);
      return;
   end
   if ic == 5
      if n1 <= th*ic*ns
         eff = n1/ic/ns;
         disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
         close(h)
         newCut = input('How many extra cut do you want to include? \n (negative value to exit) \n' );
         if newCut > 0
            Dinfo.D = A';
            Dinfo.uq = uq;
            opt.SampleOption.ExtraCut = ceil(newCut);
            [xVal2,eff] = RSP(obj,N-n1,xS,Dinfo,opt);
            xVal = [xVal; xVal2];
%             xVal = xVal(randperm(size(xVal,1),size(xVal,1)),:);
            return;
         else
            return;
         end
      end
   end
end

