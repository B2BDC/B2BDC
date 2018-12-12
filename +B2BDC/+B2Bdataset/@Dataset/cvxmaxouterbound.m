function [maxout,s,xs] = cvxmaxouterbound(obj,q,frac,abE,rflag)
% Calculate the outer bound of maximum values of the dataset
% target QOI subject to the constraints within the dataset through cvx
% QOIobj - A B2BDC.B2Bmodel.Model object specifies the QOI
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: June 15, 2015   Wenyu Li
%  Modified: July 5, 2015    Wenyu Li (Sensitivity added)
%  Modified: July 14, 2015    Wenyu Li (Extra combination for linear
%                                       constraints added)

units = obj.DatasetUnits.Values;
bds = obj.calBound;
bds(:,1) = bds(:,1) - abE;
bds(:,2) = bds(:,2) + abE;
[idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion(q);
id0 = idall{end};
Q0 = Qall{end};
N0 = Nall{end};
D0 = Dall{end};
name2 = obj.VarNames;
vu = [obj.Variables.Values.UpperBound]';
vl = [obj.Variables.Values.LowerBound]';
if obj.ModelDiscrepancyFlag
   nMD = obj.ModelDiscrepancy.Variables.Length;
   HMD = obj.ModelDiscrepancy.Variables.calBound;
   vu = [vu; HMD(:,2)];
   vl = [vl; HMD(:,1)];
else
   nMD = 0;
   HMD = [];
end
if obj.ParameterDiscrepancyFlag
   nPD = obj.ParameterDiscrepancy.Variables.Length;
   HPD = obj.ParameterDiscrepancy.Variables.calBound;
   vu = [vu; HPD(:,2)];
   vl = [vl; HPD(:,1)];
else
   nPD = 0;
   HPD = [];
end
vd = vu - vl;
if rflag
   error('Predicted QOI has variables not in the dataset');
else
   n_units = obj.Length;
end
n_variable = obj.Variables.Length;
[Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ,LBD]  = obj.getInequalQuad(bds,frac);
d = bds(:,2)-bds(:,1);
nL = length(Qx) - n_variable - nMD -nPD;
if nL > 0
   tlb = LBD(:,1);
   tub = LBD(:,2);
   tdb = tub-tlb;
   L2 = L(n_variable+nMD+nPD+1:end,:);
end
cvx_begin sdp quiet
cvx_precision('best');
variable tmax
variable lamEL(n_units)
variable lamEU(n_units)
variable lamV(n_variable+nMD+nPD)
if nL > 0
   variable lamL(nL)
end
if n_extra ~= 0
   variable lamExtra(n_extra)
end
dual variable S
if ~isempty(N0)
   targetN = zeros(n_variable+nMD+nPD+1);
   targetD = zeros(n_variable+nMD+nPD+1);
   targetN([1;id0+1],[1;id0+1]) = N0;
   targetD([1;id0+1],[1;id0+1]) = D0;
   cvxconstr = tmax*targetD-targetN;
elseif ~isempty(Q0)
   targetMatrix1 = zeros(n_variable+nMD+nPD+1);
   targetMatrix2 = zeros(n_variable+nMD+nPD+1);
   targetMatrix2(1,1) = 1;
   targetMatrix1([1;id0+1],[1;id0+1]) = Q0;
   cvxconstr = tmax*targetMatrix2-targetMatrix1;
end
for i = 1:n_units
   cvxconstr = cvxconstr + lamEU(i)*Qunits{2*i-1}+...
      lamEL(i)*Qunits{2*i};
end
for i = 1:n_variable+nMD+nPD
   cvxconstr = cvxconstr + lamV(i)*Qx{i};
end
if nL > 0
   for i = 1:nL
      cvxconstr = cvxconstr + lamL(i) * Qx{n_variable+nMD+nPD+i};
   end
end
if n_extra ~= 0
   for k = 1:n_extra
      cvxconstr = cvxconstr + lamExtra(k)*Qextra{k};
   end
end
minimize tmax
subject to
S: cvxconstr >= zeros(n_variable+nMD+nPD+1,n_variable+nMD+nPD+1);
lamEL >= 0;
lamEU >= 0;
lamV >= 0;
if nL > 0
   lamL >= 0
end
if n_extra ~= 0
   lamExtra >= 0;
end
cvx_end
maxout = cvx_optval;
xs = S(2:n_variable+nMD+nPD+1,1);
s.expu = zeros(n_units,1);
s.expl = zeros(n_units,1);
s.expu(~idRQ) = lamEU(~idRQ).*d(~idRQ)*S(1,1);
s.expl(~idRQ) = lamEL(~idRQ).*d(~idRQ)*S(1,1);
for j1 = 1:sum(idRQ)
   ida = find(idRQ,j1);
   idU = ida(end);
   D = units(idU).SurrogateModel.Denominator;
   name1 = units(idU).SurrogateModel.VarNames;
   [~,~,id] = intersect(name1,name2,'stable');
   id = [1;id+1];
   s.expu(idU) = lamEU(idU)*d(idU)*trace(S(id,id)*D');
   s.expl(idU) = lamEL(idU)*d(idU)*trace(S(id,id)*D');
end
s.varu = (-vl*S(1,1) + S(2:end,1)).*lamV.*vd;
s.varl = (vu*S(1,1) - S(2:end,1)).*lamV.*vd;
if nL > 0
   L2S = L2*S(2:end,1);
   s.linu = lamL.*(-tlb*S(1,1) + L2S).*tdb;
   s.linl = lamL.*(tub*S(1,1) - L2S).*tdb;
else
   s.linu = [];
   s.linl = [];
end
if n_extra ~= 0
   LS = L*S(2:end,1);
   for k1 = 1:n_extra
      i1 = extraIdx(k1,1);
      j1 = extraIdx(k1,1);
      c1 = extraIdx(k1,1);
      if c1 == 1
         if i1 <= n_variable
            li = vl(i1);
            if j1 <= n_variable
               lj = vl(j1);
               s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            else
               lj = tlb(j1-n_variable);
               s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            end
            s.varl(i1) = s.varl(i1) + vd(i1)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
         else
            li = tlb(i1-n_variable);
            if j1 <= n_variable
               lj = vl(j1);
               s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            else
               lj = tlb(j1-n_variable);
               s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            end
            s.linl(i1-n_variable) = s.linl(i1-n_variable) + ...
               tdb(i1-n_variable)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
         end
      elseif c1 == 2
         if i1 <= n_variable
            li = vl(i1);
            if j1 <= n_variable
               uj = vu(j1);
               s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            else
               uj = tub(j1-n_variable);
               s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            end
            s.varl(i1) = s.varl(i1) + vd(i1)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
         else
            li = tlb(i1-n_variable);
            if j1 <= n_variable
               uj = vu(j1);
               s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            else
               uj = tub(j1-n_variable);
               s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
            end
            s.linl(i1-n_variable) = s.linl(i1-n_variable) + ...
               tdb(i1-n_variable)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
         end
      elseif c1 == 3
         if i1 <= n_variable
            ui = vu(i1);
            if j1 <= n_variable
               lj = vl(j1);
               s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            else
               lj = tlb(j1-n_variable);
               s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            end
            s.varu(i1) = s.varu(i1) + vd(i1)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
         else
            ui = tub(i1-n_variable);
            if j1 <= n_variable
               lj = vl(j1);
               s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            else
               lj = tlb(j1-n_variable);
               s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            end
            s.linu(i1-n_variable) = s.linu(i1-n_variable) + ...
               tdb(i1-n_variable)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
         end
      else
         if i1 <= n_variable
            ui = vu(i1);
            if j1 <= n_variable
               uj = vu(j1);
               s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            else
               uj = tub(j1-n_variable);
               s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            end
            s.varu(i1) = s.varu(i1) + vd(i1)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
         else
            ui = tub(i1-n_variable);
            if j1 <= n_variable
               uj = vu(j1);
               s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            else
               uj = tub(j1-n_variable);
               s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
                  tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
            end
            s.linu(i1-n_variable) = s.linu(i1-n_variable) + ...
               tdb(i1-n_variable)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
         end
      end
   end
end
