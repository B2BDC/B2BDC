function newVar = calScale(obj)

% S = CALSCALE(OBJ) returns a nVar-by-1 sclae vector such that the hit and
% run algorithm works most effective if each variable is scaled as
% suggested by s.

%  Created: August 29, 2016     Wenyu Li

newVar = obj;
if isempty(obj.ScalingVector)
   tolerance = 1e-5;
   nVar = obj.Length;
   xbd = obj.calBound;
   vd = xbd(:,2)-xbd(:,1);
   Ax = diag(1./vd);
   Ax = [-Ax, ones(nVar,1), zeros(nVar,1);
      Ax, zeros(nVar,1), -ones(nVar,1)];
   bx = zeros(2*nVar,1);
   A = obj.ExtraLinConstraint.A;
   if ~isempty(A)
      LB = obj.ExtraLinConstraint.LB;
      UB = obj.ExtraLinConstraint.UB;
      eqTest = (UB-LB)./sum(abs(A),2);
      ieq = (eqTest <= tolerance);
      nA = length(LB) - sum(ieq);
      dA = UB-LB;
      if nA > 0
         AL = A(~ieq,:);
         AL = [-AL, dA(~ieq).*ones(nA,1), zeros(nA,1);
            AL, zeros(nA,1), -dA(~ieq).*ones(nA,1)];
         bL = zeros(2*nA,1);
      end
   else
      AL = [];
      bL = [];
   end
   A = [Ax;AL];
   b = [bx;bL];
   f = [zeros(nVar,1);-1;1];
   lb = [ones(nVar,1);-inf(2,1)];
   ub = inf(nVar+2,1);
   warning('off','all');
   opt = optimoptions('linprog');
   opt.Display = 'none';
   [x,feval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],opt);
   warning('on','all');
   if exitflag < 0
      s = ones(nVar,1);
   else
      s = x(1:nVar)/min(x(1:nVar));
   end
   newVar.ScalingVector = s;
end