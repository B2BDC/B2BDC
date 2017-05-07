function vars = pruneLinConstraint(obj)

% Check on each linear constraint, delete the ones that are redundant

%  Created: Oct 13, 2016   Wenyu Li

A = obj.ExtraLinConstraint.A;
LB = obj.ExtraLinConstraint.LB;
UB = obj.ExtraLinConstraint.UB;
n1 = size(A,1);
xbd = obj.calBound;
vars = obj;
vars.clearExtraConstraint;
ic = 0;
cri = 1e-4;
warning('off','all');
for i = 1:n1
   tA = A;
   tUB = UB;
   tLB = LB;
   tA(i-ic,:) = [];
   tUB(i-ic) = [];
   tLB(i-ic) = [];
   opt = optimoptions('linprog');
   opt.Display = 'none';
   [~,lb] = linprog(A(i-ic,:)',[tA;-tA],[tUB;-tLB],[],[],...
      xbd(:,1),xbd(:,2),[],opt);
   [~,ub] = linprog(-A(i-ic,:)',[tA;-tA],[tUB;-tLB],[],[],...
      xbd(:,1),xbd(:,2),[],opt);
   ub = -ub;
   if lb >= LB(i-ic)-cri && ub <= UB(i-ic)+cri
      A(i-ic,:) = [];
      LB(i-ic) = [];
      UB(i-ic) = [];
      ic = ic+1;
   end
end
warning('on','all');
vars.addLinearConstraint([A;-A],[UB;-LB]);