function [xVal,eff] = RSH(obj,N,xF,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 3, 2016     Wenyu Li

xVal = [];
sOpt = opt.SampleOption;
th2 = sOpt.RejectionTol;
s1 = sOpt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nPC = nVar - sOpt.TruncatedPC;
if nPC < 1
   nPC = 1;
end
% t1 = toc;
if isempty(Dinfo)
   if isempty(xF)
      nS = min(10^7,10*nVar^3);
      xF = obj.collectSamples(nS,[],opt);
   else
      nS = size(xF,1);
   end
end
if isempty(Dinfo)
   xAve = mean(xF);
   xC = xF - repmat(xAve,nS,1);
   [V,D] = eig((xC' * xC)/nS);
   Ddiag = diag(D);
   [~,id] = sort(Ddiag,'descend');
   vv = V(:,id);
else
   uqD = Dinfo.uq;
   xAve = 0.5*(uqD(:,1)+uqD(:,2))';
   vv = Dinfo.D;
   vv = normc(vv);
   dd = uqD(:,2) - uqD(:,1);
   [~,id] = sort(dd,'descend');
   vv = vv(:,id);
   uqD = uqD - [xAve' xAve'];
   uqD = uqD(id,:);
end
if isempty(Dinfo)
   uq = zeros(nPC,2);
   switch s1
      case 'Outer'
         opt.Prediction = 'outer';
         opt.ExtraLinFraction = -1;
         %       h = waitbar(0,'Calculating directional uncertainty...');
         for i = 1:nPC
            tv = [-xAve*vv(:,i); vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uq(i,1) = preQ.min;
            uq(i,2) = preQ.max;
            %          waitbar(i/nPC,h);
         end
         %       close(h);
      case 'Inner'
         opt.Prediction = 'inner';
         %       h = waitbar(0,'Calculating directional uncertainty...');
         for i = 1:nPC
            tv = [-xAve*vv(:,i); vv(:,i)];
            tM = generateModel(tv,vars);
            preQ = obj.predictQOI(tM,opt);
            uq(i,1) = preQ.min;
            uq(i,2) = preQ.max;
            %          waitbar(i/nPC,h);
         end
         %       close(h);
      case 'Sample'
         %       h = waitbar(0,'Calculating directional uncertainty...');
         for i = 1:nPC
            f = xC*vv(:,i);
            uq(i,1) = min(f);
            uq(i,2) = max(f);
            %          waitbar(i/nPC,h);
         end
         %       close(h);
      case 'Truncation'
         %       h = waitbar(0,'Calculating directional uncertainty...');
         sf = 1-sOpt.UncertaintyTruncation;
         for i = 1:nPC
            f = xC*vv(:,i);
            uq(i,1) = sf*min(f);
            uq(i,2) = sf*max(f);
            %          waitbar(i/nPC,h);
         end
         %       close(h);
   end
else
   uq = uqD(1:nPC,:);
end
% t2 = toc;
% tt = t2-t1;
ic = 1;
n1 = size(xVal,1);
ns = sOpt.BatchMaxSample;
% h = waitbar(n1/N,'Collecting uniform samples of the feasible set...');
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
%    waitbar(min(n1/N,1),h);
   if n1 >= N
      eff = n1/ic/ns;
%       close(h);
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
xVal = xVal(1:N,:);