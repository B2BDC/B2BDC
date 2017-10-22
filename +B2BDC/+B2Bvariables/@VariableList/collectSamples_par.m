function Xvals = collectSamples_par(obj,N,xStart,opt)
%   XVALS = COLLECTSAMPLES(OBJ, N, XSTART,OPT) returns a n-by-nVariable matrix of
%   feasible points of the domain specified by OBJ, where N is the number of
%   feasible points to be returned.

if nargin < 4
   opt = generateSampleOpt;
end
cMax = 10^4;
tolerance = 1e-5;
if opt.ParameterScaling
   obj = obj.calScale;
   s = obj.ScalingVector;
else
   s = ones(obj.Length,1);
end
A0 = obj.ExtraLinConstraint.A;
Q = obj.ExtraQuaConstraint.Q;
Qb = obj.ExtraQuaConstraint.UB;
if isempty(Q)
   if ~isempty(A0)
      UB = obj.ExtraLinConstraint.UB;
      LB = obj.ExtraLinConstraint.LB;
      eqTest = (UB-LB)./sum(abs(A0),2);
      A0 = A0.*repmat(s',size(A0,1),1);
      ieq = (eqTest <= tolerance);
      if any(ieq)
         Neq = sum(ieq);
         A1 = [A0(~ieq,:); -A0(~ieq,:)];
         B1 = [UB(~ieq); -LB(~ieq)];
         Aeq = A0(ieq,:);
         Beq = 0.5*(LB(ieq)+UB(ieq));
      else
         A1 = [A0; -A0];
         Aeq = [];
         B1 = [UB; - LB];
         Beq = [];
      end
   else
      A1 = [];
      B1 = [];
      Aeq = [];
      Beq = [];
   end
   xbd = obj.calBound;
   xLB = xbd(:,1);
   xUB = xbd(:,2);
   A2 = [diag(s); -diag(s)];
   B2 = [xUB; -xLB];
   Arw = [A1; A2];
   brw = [B1; B2];
   if ~isempty(Aeq)
      bs1 = orth(Aeq');
      bs2 = null(Aeq);
      Tm = [bs2, bs1]';
   end
   
   if nargin < 3 || isempty(xStart)
      iN = (s<0);
      tLB = xLB./s;
      tUB = xUB./s;
      tt = tLB(iN);
      tLB(iN) = tUB(iN);
      tUB(iN) = tt;
      if isempty(A0)
         x0 = obj.makeLHSsample(1);
         x0 = x0./s;
         cc = 0;
         while ~all(Arw*x0 < brw) && cc < cMax
            disp('LOL')
            x0 = obj.makeLHSsample(1);
            x0 = x0./s;
            cc = cc+1;
         end
         if cc == cMax
            x0 = [];
         end
      else
         warning('off','all');
         opt1 = optimoptions('linprog');
         opt1.Display = 'none';
         x0 = linprog(zeros(size(Arw,2),1),...
            Arw,brw,Aeq,Beq,tLB,tUB,opt1);
         cc = 0;
         while ~all(Arw*x0 < brw) && cc < cMax
            disp('LOL')
            x0 = linprog(zeros(size(Arw,2),1),...
               Arw,brw,Aeq,Beq,tLB,tUB,opt1);
         end
         if cc == cMax
            x0 = [];
         end
      end
   else
      if obj.isFeasiblePoint(xStart)
         x0 = xStart./s;
      else
         iN = s<0;
         tLB = xLB./s;
         tUB = xUB./s;
         tt = tLB(iN);
         tLB(iN) = tUB(iN);
         tUB(iN) = tt;
         warning('off','all');
         opt1 = optimoptions('linprog');
         opt1.Display = 'none';
         x0 = linprog(zeros(size(Arw,2),1),...
            Arw,brw,Aeq,Beq,tLB,tUB,opt1);
      end
   end
   warning('on','all');
   
   if ~isempty(x0)
      xInit = x0;
   else
      errordlg('Infeasible problem: the domain defined cannot be shown to be feasible');
      error('collectSamples:Inconclusive',...
         'Domain cannot be shown to be Consistent');
   end
   
   if ~isempty(Aeq)
      y0 = Tm*x0;
      Ay = Arw*Tm';
      by = brw - Ay(:,end-Neq+1:end) * y0(end-Neq+1:end);
   end
   % Check Dataset model type
   if isempty(Aeq)
      Xvals = HR(Arw, brw, xInit);
   else
      Yvals = HR(Ay(:,1:end-Neq), by, y0(1:end-Neq));
      Yvals = [Yvals; repmat(y0(end-Neq+1:end),1,N)];
      Xvals = Tm' * Yvals;
   end
   Xvals = Xvals';
   Xvals = Xvals .* repmat(s',N,1);
else
   H = obj.calBound;
   nVar = size(H,1);
   nQ = length(Q);
   nL = size(A0,1);
   UB = obj.ExtraLinConstraint.UB;
   LB = obj.ExtraLinConstraint.LB;
   if nargin >= 3 && ~isempty(xStart) && obj.isFeasiblePoint(xStart)
      xInit = xStart;
   else
      xInit = obj.ExtraQuaConstraint.xStart;
   end
   Mgd = cell(nQ+nL,1);
   idx = cell(nQ+nL,1);
   for i1 = 1:nQ
      idx{i1} = (1:nVar)';
      Mgd{i1} = Q{i1};
      Mgd{i1}(1,1) = Mgd{i1}(1,1) - Qb(i1);
   end
   for i1 = 1:nL
      ai = A0(i1,:);
      ub = UB(i1);
      lb = LB(i1);
      tmpM = zeros(nVar+1);
      tmpM(2:end,2:end) = ai' * ai;
      tmpM(1) = ub*lb;
      tmpM(1,2:end) = -0.5*(ub+lb)*ai;
      tmpM(2:end,1) = -0.5*(ub+lb)*ai';
      Mgd{i1+nQ} = tmpM;
      idx{i1+nQ} = (1:nVar)';
   end
   Xvals = zeros(N, nVar);
   ic = 1;
   nstep = opt.StepInterval;
   n2 = nstep*nVar;
   for i1=1:n2
      V = randn(1, nVar);
      zPDir = find(xInit >= H*[0.01; 0.99]);
      zNDir = find(xInit <= H*[0.99; 0.01]);
      V(zPDir) = -abs(V(zPDir));
      V(zNDir) = abs(V(zNDir));
      tmpx = B2BDC.B2Bdataset.Dataset.q2sample(Mgd,idx,H,xInit,V');
      if mod(i1,nstep) == 0
         Xvals(ic,:) = tmpx';
         ic = ic+1;
      end
      xInit = tmpx;
   end
end


   function X = HR(A,b,xFeas)
      % Hit-and-Run for polytope described by Ax<=b, where A is m-by-n, and b is
      % m-by 1.   Code assumes polytope is bounded.  xFeas is n-by-1, and is
      % strictly feasible point.   N is the desired number of samples.  The
      % output, X, is n-by-N
      
      n1 = size(A,2);
      if ~all(A*xFeas <= b)
         count = 0;
         flag = 0;
         while count < 100
            xt = xFeas + 1e-10 * randn(n1,1);
            count = count+1;
            if all(A*xt <= b)
               xFeas = xt;
               flag = 1;
               break
            end
         end
         if flag == 0
            error('The starting point is not feasible')
         end
      end
      nstep = opt.StepInterval;
      n2 = nstep*(N-1)+1;
      nss = round(10^7/size(A,1));
      ns1 = floor(n2/nss)+1;
      ns2 = mod(n2,nss);
      X = zeros(n1,N);
      ic = 1;
      it = 1;
      for k = 1:ns1
         if k == ns1
            n2 = ns2;
         else
            n2 = nss;
         end
         V = randn(n1,n2);
         R = rand(1,n2);
         AV = A*V;
         if k == 1
            xt = xFeas;
         end
         for i=1:n2
            bmAx = b - A*xt;
            AVi = AV(:,i);
            pos = AVi>0;
            neg = AVi<0;
            tValues = bmAx./AVi;
            tmax = min(tValues(pos));
            tmin = max(tValues(neg));
            t = R(i)*tmin + (1-R(i))*tmax;
            if mod(it,nstep) == 1
               X(:,ic) = xt + t*V(:,i);
               ic = ic+1;
            end
            xt = xt + t*V(:,i);
            it = it+1;
         end
      end
      X = X(:,randperm(N,N));
   end
end