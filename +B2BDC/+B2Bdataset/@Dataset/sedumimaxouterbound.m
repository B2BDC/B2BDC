function [maxout,s] = sedumimaxouterbound(obj,q, frac, abE, rflag)
% Calculate the outer bound of maximum values of the dataset
% target QOI subject to the constraints within the dataset through sedumi
% QOIobj - A B2BDC.B2Bmodel.Model object specifies the QOI
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015   Wenyu Li
%  Modified: Feb 16, 2016  Wenyu Li (Normalized sensitivity calculated) 

bds = obj.calBound;
bds(:,1) = bds(:,1) - abE;
bds(:,2) = bds(:,2) + abE;
if rflag
   error('Predicted QOI has variables not in the dataset');
else
   n_units = obj.Length;
end
units = obj.DatasetUnits.Values;
name2 = obj.VarNames;
[idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion(q);
id0 = idall{end};
Q0 = Qall{end};
N0 = Nall{end};
D0 = Dall{end};
[Qunits, Qx, Qextra, n_extra, extraIdx, L, idRQ, LBD]  = obj.getInequalQuad(bds,frac);
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
n_variable = obj.Variables.Length;
[Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
bds = bds(1:n_units,:);
d = bds(:,2)-bds(:,1);
s = [];
nL = length(Qx) - n_variable -nMD-nPD;
n_opt = 1+2*n_units+n_variable+nMD+nPD+nL+n_extra;
if nL > 0
   tlb = LBD(:,1);
   tub = LBD(:,2);
   tdb = tub-tlb;
   L2 = L(n_variable+nMD+nPD+1:end,:);
end

[A,b,c,K] = setSedumi;
pars.fid = 0;
[xopt,yopt,info] = sedumi(A,b,c,K,pars);
S = reshape(xopt(n_opt:end),n_variable+nMD+nPD+1,n_variable+nMD+nPD+1);
S = 0.5*(S+S');

if info.pinf ~= 0 || info.dinf ~= 0
   disp('Not both primal/dual problem are feasible, please review the results with cautions')
end
maxout = yopt(end);
lamEU = yopt(1:n_units);
lamEL = yopt(n_units+1:2*n_units);
lamV = yopt(2*n_units+1:2*n_units+n_variable+nMD+nPD);
if nL > 0
   lamL = yopt(2*n_units+n_variable+nMD+nPD+1:2*n_units+n_variable+nMD+nPD+nL);
end
if n_extra ~= 0
   lamExtra = yopt(2*n_units+nL+n_variable+nMD+nPD+1:end-1);
end
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




   function [A,b,c,K] = setSedumi()
      Alam = [-speye(n_opt-1), spalloc(n_opt-1,1,0)];
      clam = spalloc(n_opt-1,1,0);
      As = zeros((n_variable+nMD+nPD+1)^2,n_opt);
      if ~isempty(N0)
         targetN = zeros(n_variable+nMD+nPD+1);
         targetD = zeros(n_variable+nMD+nPD+1);
         targetN([1;id0+1],[1;id0+1]) = N0;
         targetD([1;id0+1],[1;id0+1]) = D0;
         targetMatrix = -targetN;
      elseif ~isempty(Q0)
         targetMatrix = zeros(n_variable+nMD+nPD+1);
         targetMatrix([1;id0+1],[1;id0+1]) = -Q0;
      end
      cs = targetMatrix(:);
      for i = 1:n_units
         As(:,i) = -Qunits{2*i-1}(:);
         As(:,i+n_units) = -Qunits{2*i}(:);
      end
      for i = 1:n_variable+nMD+nPD
         As(:,2*n_units+i) = -Qx{i}(:);
      end
      for i = 1:nL
         As(:,2*n_units+n_variable+nMD+nPD+i) = -Qx{n_variable+nMD+nPD+i}(:);
      end
      if n_extra ~= 0
         for k = 1:n_extra
            As(:,2*n_units+n_variable+nMD+nPD+nL+k) = -Qextra{k}(:);
         end
      end          
      if ~isempty(N0)
         As(:,end) = -targetD;
      elseif ~isempty(Q0)
         As(1,end) = -1;
      end
      At = [Alam; As];
      c = [clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = -1;
      K.l = n_opt-1;
      K.s = n_variable+nMD+nPD+1;
      A = At';    
   end
end