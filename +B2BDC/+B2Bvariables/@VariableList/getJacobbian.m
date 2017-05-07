function J = getJacobbian(obj,J0)
% J=GETJACOBBIAN(OBJ,J0) returns a information matrix for the selection of
% extra linear constraint combinations used in optimization

%  Created: April 29, 2016     Wenyu Li

A = obj.ExtraLinConstraint.A;
LB = abs(obj.ExtraLinConstraint.LB);
UB = abs(obj.ExtraLinConstraint.UB);
nA = size(A,1);
J = zeros(nA);
for i = 1:nA-1
   for j = i+1:nA
      ai = abs(A(i,:));
      aj = abs(A(j,:));
      c = 0.25*(1/LB(i)/LB(j) + 1/LB(i)/UB(j) + 1/UB(i)/LB(j) + 1/UB(i)/UB(j));
      J(i,j) = c*trace(ai'*aj * J0);
   end
end