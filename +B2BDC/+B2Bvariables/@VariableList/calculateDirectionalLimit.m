function dRange = calculateDirectionalLimit(obj,d)

% DRANGE = CALCULATEDIRECTIONALLIMIT(OBJ,D) calculates the directional limit
% along the direction d of the domain specified by obj. It assumes the
% direction passing through the origin, specified by the nominal values of
% obj. The direction vector is normalized.

%  Created: Jan 17, 2018   Wenyu Li

nd = size(d,2);
d = normc(d);
dRange = zeros(nd,2);
H = obj.calBound;
x0 = obj.calNominal;
if ~isempty(obj.ExtraLinConstraint.A)
   A = obj.ExtraLinConstraint.A;
   LB = obj.ExtraLinConstraint.LB;
   UB = obj.ExtraLinConstraint.UB;
   dH1 = [H(:,2)-x0; UB-A*x0];
   dH2 = [x0-H(:,1); A*x0-LB];
else
   dH1 = H(:,2)-x0;
   dH2 = x0-H(:,1);
end
for i = 1:nd
   dd = d(:,i);
   if ~isempty(obj.ExtraLinConstraint.A)
      dd = [dd;A*dd];
   end
   n1 = length(dd);
   idPos = dd>=0;
   idNeg = dd<0;
   tt = inf(n1,1);
   tt(idPos) = dH2(idPos)./dd(idPos);
   tt(idNeg) = -dH1(idNeg)./dd(idNeg);
   dRange(i,2) = min(tt);
   tt = -inf(n1,1);
   tt(idPos) = -dH1(idPos)./dd(idPos);
   tt(idNeg) = dH2(idNeg)./dd(idNeg);
   dRange(i,1) = max(tt);
end