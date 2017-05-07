function sparsePoly = createSparsePOP(obj,flag)
% SPARSEPOLY = CREATESPARSEPOP(OBJ,FLAG) creates a sparsePOP format object from
% the PolyModel object. The sparsePOP object is compatible with SparsePOP
% toolbox. The input flag is used to denote the object type: 1 for
% objective or inequality constraint and -1 for equality constraint.

%  Created: Dec 15, 2016     Wenyu Li

if nargin > 1
   sparsePoly.typeCone = flag;
else
   sparsePoly.typeCone = 1;
end
s1 = obj.SupportMatrix;
sparsePoly.dimVar = size(s1,2);
sparsePoly.degree = obj.Degree;
sparsePoly.noTerms = size(s1,1);
sparsePoly.supports = s1;
sparsePoly.coef = obj.Coefficient;
