function y = isFeasiblePoint(obj, x)
% Y=ISFEASIBLEPOINT(OBJ,X) returns an array of logical results specifying
% whether the corresponding point is feasible. The input x should be a
% nSample-by-nVar size matrix.

%  Created: June 28, 2016     Wenyu Li

tolerance = 1e-12;
if size(x,2) ~= obj.Length
   x = x';
end
if size(x,2) ~= obj.Length
   error('Invalid input sample point size')
end
nS = size(x,1);
y = false(nS,1);
H = obj.calBound;
xlb = H(:,1);
xub = H(:,2);
A = obj.ExtraLinConstraint.A;
LB = obj.ExtraLinConstraint.LB;
UB = obj.ExtraLinConstraint.UB;
res1 = repmat(xlb,1,nS) <= x' & x'<= repmat(xub,1,nS);
if ~isempty(A)
   res2 = repmat(LB-tolerance,1,nS) <= A*x' & repmat(UB+tolerance,1,nS) >= A*x';
   res = [res1; res2];
else
   res = res1;
end
id = all(res);
y(id) = true;