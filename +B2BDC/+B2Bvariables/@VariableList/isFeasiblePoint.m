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
H = obj.calBound;
xlb = H(:,1);
xub = H(:,2);
A = obj.ExtraLinConstraint.A;
LB = obj.ExtraLinConstraint.LB;
UB = obj.ExtraLinConstraint.UB;
r1 = true(nS,1);
n1 = ceil(nS/1e7);
d1 = mod(nS,1e7);
if d1 == 0
   d1 = 1e7;
end
count = 0;
for i = 1:n1
   if i < n1
      t1 = all(repmat(xlb,1,1e7) <= x(count+1:count+1e7,:)',1);
      t2 = all(x(count+1:count+1e7,:)'<= repmat(xub,1,1e7),1);
      r1(count+1:count+1e7) = all([t1;t2],1);
      count = count+1e7;
   else
      t1 = all(repmat(xlb,1,d1) <= x(count+1:count+d1,:)',1);
      t2 = all(x(count+1:count+d1,:)'<= repmat(xub,1,d1),1);
      r1(count+1:end) = all([t1;t2],1);
   end
end
% res1 = repmat(xlb,1,nS) <= x' & x'<= repmat(xub,1,nS);
if ~isempty(A)
   r2 = true(nS,1);
   n22 = round(1e8/size(A,1));
   n2 = ceil(nS/n22);
   d2 = mod(nS,n22);
   if d2 == 0
      d2 = n22;
   end
   count = 0;
   for i = 1:n2
      if i < n2
         t1 = all(repmat(LB-tolerance,1,n22) <= A*x(count+1:count+n22,:)',1);
         t2 = all(A*x(count+1:count+n22,:)'<= repmat(UB+tolerance,1,n22),1);
         r2(count+1:count+n22) = all([t1;t2],1);
         count = count+n22;
      else
         t1 = all(repmat(LB-tolerance,1,d2) <= A*x(count+1:count+d2,:)',1);
         t2 = all(A*x(count+1:count+d2,:)'<= repmat(UB+tolerance,1,d2),1);
         r2(count+1:end) = all([t1;t2],1);
      end
   end
   y = all([r1,r2],2);
else
   y = r1;
end