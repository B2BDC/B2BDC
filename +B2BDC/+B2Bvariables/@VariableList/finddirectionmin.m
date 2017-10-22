function y = finddirectionmin(obj,A,b,c)
% Y = FINDDIRECTIONMIN(OBJ,A,B,C) returns the redundant maxmum lower bound
% corresponds to the linear constraint c'*x in the VariableList
% object.

%  Created: May 3, 2016     Wenyu Li

flag = ~isempty(obj.ExtraLinConstraint.A);
if flag
   Ai = [obj.ExtraLinConstraint.A; -obj.ExtraLinConstraint.A];
   bi = [obj.ExtraLinConstraint.UB; -obj.ExtraLinConstraint.LB];
   ni = length(bi);
end
H = obj.calBound;
nVar = obj.Length;
n1 = size(A,1);
if flag
   nS = n1 + ni + 2*nVar;
else
   nS = n1 + 2*nVar;
end
A0 = zeros(nS,nVar+nS);
A0(1:n1,1:nVar) = A;
if flag
   A0(n1+1:n1+ni,1:nVar) = Ai;
end
A0(1:nS,nVar+1:end) = eye(nS);
if flag
   A0(n1+ni+1:end,1:nVar) = [eye(nVar); -eye(nVar)];
else
   A0(n1+1:end,1:nVar) = [eye(nVar); -eye(nVar)];
end
if flag
   b = [b; bi; H(:,2); -H(:,1)];
else
   b = [b; H(:,2); -H(:,1)];
end
K.f = nVar;
K.l = nS;
cc = zeros(nVar+nS,1);
cc(1:nVar) = c';
pars.fid = 0;
[x,~,~] = sedumi(A0,b,cc,K,pars);
y = c*x(1:nVar);

