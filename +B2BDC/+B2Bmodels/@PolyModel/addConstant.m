function newPoly = addConstant(obj,c)
% NEWPOLY = ADDCONSTANT(OBJ,C) returns a new PolyModel object by adding a
% constant term c to the original PolyModel obj.

%  Created: Dec 16, 2016     Wenyu Li

newPoly = obj;
s1 = newPoly.SupportMatrix;
[~,id] = intersect(s1,zeros(1,size(s1,2)),'rows');
if ~isempty(id)
   newPoly.Coefficient(id) = newPoly.Coefficient(id)+c;
else
   s1 = [s1; zeros(1,size(s1,2))];
   newPoly.SupportMatrix = s1;
   newPoly.Coefficient(end+1) = c;
end