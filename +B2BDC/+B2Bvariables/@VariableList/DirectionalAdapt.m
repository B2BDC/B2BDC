function [sVar,flag] = DirectionalAdapt(obj,x,ss,monoc)

% [SVAR,FLAG] = DIRECTIONALADAPT(OBJ,X,SS) updates the polytope defined by
% obj. The update is based on the provided samples x, so the regions of
% negligable directional density are cut while boundaries with observable
% density are expanded.

%  Created: April 23, 2018  Wenyu Li

if nargin<4
   monoc = true;
end
if nargin<3 || ss<0 || ss>1
   ss = 0;
end
nS = size(x,1);
nbin = max(round(0.02*nS),20);
A = obj.ExtraLinConstraint.A;
LB = obj.ExtraLinConstraint.LB;
UB = obj.ExtraLinConstraint.UB;
nD = size(A,1);
rr = x*A';
dR = (max(rr)-min(rr))/nbin;
tLB = LB+dR';
tUB = UB-dR';
tVar = obj.clearExtraConstraint;
tVar = tVar.addLinearConstraint(A,tLB,tUB);
tS = sum(tVar.isFeasiblePoint(x))/nS;
flag = false(nD,1);
if tS <= ss
   sVar = obj;
   return
end
ss = ss/nD;
nL = sum(rr<=repmat(LB'+dR,nS,1))/nS;
nR = sum(rr>=repmat(UB'-dR,nS,1))/nS;
if monoc
   for i = 1:nD
      if nL(i) > ss
         LB(i) = LB(i)-2.5*dR(i);
         flag(i) = true;
      else
         LB(i) = min(LB(i),quantile(rr(:,i),0.5*ss)-1.5*dR(i));
      end
      if nR(i) > ss
         UB(i) = UB(i)+2.5*dR(i);
         flag(i) = true;
      else
         UB(i) = max(UB(i),quantile(rr(:,i),1-0.5*ss)+1.5*dR(i));
      end
   end
else
   for i = 1:nD
      if nL(i) > ss
         LB(i) = LB(i)-2.5*dR(i);
         flag(i) = true;
      else
         LB(i) = quantile(rr(:,i),0.5*ss)-1.5*dR(i);
      end
      if nR(i) > ss
         UB(i) = UB(i)+2.5*dR(i);
         flag(i) = true;
      else
         UB(i) = quantile(rr(:,i),1-0.5*ss)+1.5*dR(i);
      end
   end
end
sVar = obj.clearExtraConstraint;
sVar = sVar.addLinearConstraint(A,LB,UB);