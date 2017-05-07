function newVar = parseLinConstraint(obj)
% NEWVAR = PARSELINCONSTRAINT(OBJ) returns a new VariableList object with
% its extra linear constraint parsed such that any redundant linear
% constraints are removed.

% Created: July 14, 2016     Wenyu Li

A = obj.ExtraLinConstraint.A;
[nC,nVar] = size(A);
if nC ~= 0
   LB = obj.ExtraLinConstraint.LB;
   UB = obj.ExtraLinConstraint.UB;
   nS = 2*nVar + 2*nC;
   K.f = nVar;
   K.l = nS-2;
   flag = false(nC,1);
   pars.fid = 0;
   xbd = obj.calBound;
   A0 = [eye(nVar); -eye(nVar)];
   b0 = [xbd(:,2); -xbd(:,1)];
   A0 = [A;-A;A0];
   A0 = [A0, eye(nS)];
   b0 = [UB;-LB;b0];
   for i = 1:nC
      tmpA = A0;
      tmpA([i,i+nC],:) = [];
      tmpA(:,[i+nVar,i+nC+nVar]) = [];
      tmpb = b0;
      tmpb([i,i+nC]) = [];
      c1 = [A(i,:)';zeros(nS-2,1)];
      xopt = sedumi(tmpA,tmpb,c1,K,pars);
      if c1' * xopt < LB(i)
         break
      else
         xopt = sedumi(tmpA,tmpb,-c1,K,pars);
         if -c1' * xopt > UB(i)
            break
         else
            flag(i) = true;
         end
      end
   end
end
newVar = obj;
newVar.clearExtraConstraint;
A(flag,:) = [];
LB(flag) = [];
UB(flag) = [];
newVar.ExtraLinConstraint.A = A;
newVar.ExtraLinConstraint.UB = UB;
newVar.ExtraLinConstraint.LB = LB;