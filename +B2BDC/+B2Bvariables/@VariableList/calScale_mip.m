function newVar = calScale_mip(obj)

% S = CALSCALE(OBJ) returns a nVar-by-1 sclae vector such that the hit and
% run algorithm works most effective if each variable is scaled as
% suggested by s.

%  Created: Nov 8, 2016     Wenyu Li

newVar = obj;
if isempty(obj.ScalingVector)
   M = 10^4;
   tolerance = 1e-5;
   nVar = obj.Length;
   A = obj.ExtraLinConstraint.A;
   xbd = obj.calBound;
   vd = xbd(:,2)-xbd(:,1);
   Ax = diag(1./vd);
   if ~isempty(A)
      LB = obj.ExtraLinConstraint.LB;
      UB = obj.ExtraLinConstraint.UB;
      eqTest = (UB-LB)./sum(abs(A),2);
      ieq = (eqTest <= tolerance);
      nA = length(LB) - sum(ieq);
      dA = UB-LB;
      if nA > 0
         AL = A(~ieq,:);
         tdA = dA(~ieq);
         AL = AL./repmat(tdA,1,nVar);
      end
   else
      AL = [];
   end
   A = [Ax;AL];
   n1 = size(A,1);
   nx = 2+2*nVar+n1;
   Ai = [zeros(2*n1,1), -ones(2*n1,1), [A;-A], zeros(2*n1,n1+nVar)];
   Ai = [Ai; [ones(n1,1), zeros(n1,1), -A, -M*eye(n1), zeros(n1,nVar)]];
   Ai = [Ai; [ones(n1,1), zeros(n1,1), A, M*eye(n1,n1), zeros(n1,nVar)]];
   Ai = [Ai; [zeros(2*nVar,2), [-eye(nVar);eye(nVar)], zeros(2*nVar,n1), M*[-eye(nVar);eye(nVar)]]];
   bi = [zeros(3*n1,1); M*ones(n1,1); -ones(nVar,1); (M-1)*ones(nVar,1)];
   f = [-1;1;zeros(nx-2,1)];
   lb = [zeros(2,1);-inf(nVar,1);-0.5*ones(n1+nVar,1)];
   ub = [inf(nVar+2,1);1.5*ones(n1+nVar,1)];
   intcon = nVar+3:nx;
   warning('off','all');
   opt = optimoptions('intlinprog');
   opt.IntegerTolerance = 1e-6;
   opt.Display = 'none';
   [x,feval,exitflag] = intlinprog(f,intcon,Ai,bi,[],[],lb,ub,opt);
   warning('on','all');
   if exitflag < 0
      s = ones(nVar,1);
   else
      s = x(3:nVar+2);
   end
   newVar.ScalingVector = s;
end