<<<<<<< Updated upstream
function [minout,s] = sedumiminouterbound(obj, QOIobj, frac,abE)
% Calculate the outer bound of minimum values of the dataset
% target QOI subject to the constraints within the dataset through sedumi
% QOIobj - A B2BDC.B2Bmodel.Model object specifies the QOI
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015    Wenyu Li
%  Modified: Feb 16, 2016  Wenyu Li (Normalized sensitivity calculated)

bds = obj.calBound;
bds(:,1) = bds(:,1) - abE;
bds(:,2) = bds(:,2) + abE;
d = bds(:,2)-bds(:,1);
n_units = obj.DatasetUnits.Length;
units = obj.DatasetUnits.Values;
for n1 = 1:n_units
   my = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.dy;
   bds(n1,:) = (bds(n1,:) - my) / dy;
   d(n1) = d(n1)/dy;
end
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
n_variable = obj.Variables.Length;
vd = 2*ones(n_variable,1);
[Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
s = [];
nL = 0.5*length(Qx) - n_variable;
n_opt = 1+2*(n_units+n_variable+nL+2*n_extra);
if nL > 0
   tlb = obj.Variables.ExtraLinConstraint.LB;
   tub = obj.Variables.ExtraLinConstraint.UB;
   tdb = tub-tlb;
   L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
end

[A,b,c,K] = setSedumi;
pars.fid = 0;
[xopt,yopt,info] = sedumi(A,b,c,K,pars);
S = reshape(xopt(n_opt:end),n_variable+1,n_variable+1);
S = 0.5*(S+S');

if info.pinf ~= 0 || info.dinf ~= 0
   disp('Not both primal/dual problem are feasible, please review the results with cautions')
end
minout = yopt(end);
lamEU = -yopt(1:n_units);
lamEL = -yopt(n_units+1:2*n_units);
lamVU = -yopt(2*n_units+1:2*n_units+n_variable);
lamVL = -yopt(2*n_units+n_variable+1:2*(n_units+n_variable));
if nL > 0
   lamLU = -yopt(2*(n_units+n_variable)+1:2*(n_units+n_variable)+nL);
   lamLL = -yopt(2*(n_units+n_variable)+nL+1:2*(n_units+n_variable+nL));
end
if n_extra ~= 0
   lamExtra = -yopt(2*(n_units+nL+n_variable)+1:end-1);
end
s.expu = zeros(n_units,1);
s.expl = zeros(n_units,1);
s.expu(~idRQ) = lamEU(~idRQ).*d(~idRQ)*S(1,1);
s.expl(~idRQ) = lamEL(~idRQ).*d(~idRQ)*S(1,1);
for j1 = 1:sum(idRQ)
   ida = find(idRQ,j1);
   idU = ida(end);
   D = units(idU).SurrogateModel.NormalizedDenominator;
   s.expu(idU) = lamEU(idU)*d(idU)*dot(S,D);
   s.expl(idU) = lamEL(idU)*d(idU)*dot(S,D);
end
s.varu = (1.1*S(1,1) + S(2:end,1)).*lamVU.*vd;
s.varl = (1.1*S(1,1) - S(2:end,1)).*lamVL.*vd;
if nL > 0
   L2S = L2*S(2:end,1);
   s.linu = lamLU.*((-tlb-0.05*tdb)*S(1,1) + L2S).*tdb;
   s.linl = lamLL.*((tub+0.05*tdb)*S(1,1) - L2S).*tdb;
else
   s.linu = [];
   s.linl = [];
end
if n_extra ~= 0
   count = 1;
   for k1 = 1:n_extra
      i1 = extraIdx(k1,1);
      j1 = extraIdx(k1,2);
      if i1 <= n_variable
         if j1 <= n_variable
            s.varu(i1) = s.varu(i1) + 0.5*lamExtra(count+3)*(1-S(1,j1));
            s.varu(i1) = s.varu(i1) + 0.5*lamExtra(count+2)*(1+S(1,j1));
            s.varl(i1) = s.varl(i1) + 0.5*lamExtra(count)*(1+S(1,j1));
            s.varl(i1) = s.varl(i1) + 0.5*lamExtra(count+1)*(1-S(1,j1));
            s.varu(j1) = s.varu(j1) + 0.5*lamExtra(count+3)*(1-S(1,i1));
            s.varu(j1) = s.varu(j1) + 0.5*lamExtra(count+2)*(1+S(1,i1));
            s.varl(j1) = s.varl(j1) + 0.5*lamExtra(count)*(1+S(1,i1));
            s.varl(j1) = s.varl(j1) + 0.5*lamExtra(count+1)*(1-S(1,i1));
         else
            j1 = j1 - n_variable;
            s.varu(i1) = s.varu(i1) + lamExtra(count+3)*(tub(j1)-L2S(j1))/tdb(j1);
            s.varu(i1) = s.varu(i1) + lamExtra(count+2)*(-tlb(j1)+L2S(j1))/tdb(j1);
            s.varl(i1) = s.varl(i1) + lamExtra(count)*(-tlb(j1)+L2S(j1)/tdb(j1));
            s.varl(i1) = s.varl(i1) + lamExtra(count+1)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(j1) = s.linu(j1) + 0.5*lamExtra(count+3)*(1-S(1,i1));
            s.linu(j1) = s.linu(j1) + 0.5*lamExtra(count+2)*(1+S(1,i1));
            s.linl(j1) = s.linl(j1) + 0.5*lamExtra(count)*(1+S(1,i1));
            s.linl(j1) = s.linl(j1) + 0.5*lamExtra(count+1)*(1-S(1,i1));
         end
      else
         i1 = i1 - n_variable;
         if j1 <= n_variable
            s.linu(i1) = s.linu(i1) + 0.5*lamExtra(count+3)*(1-S(1,j1));
            s.linu(i1) = s.linu(i1) + 0.5*lamExtra(count+2)*(1+S(1,j1));
            s.linl(i1) = s.linl(i1) + 0.5*lamExtra(count)*(1+S(1,j1));
            s.linl(i1) = s.linl(i1) + 0.5*lamExtra(count+1)*(1-S(1,j1));
            s.varu(j1) = s.varu(j1) + lamExtra(count+3)*(tub(i1)-L2S(i1))/tdb(i1);
            s.varu(j1) = s.varu(j1) + lamExtra(count+2)*(-tlb(i1)+L2S(i1))/tdb(i1);
            s.varl(j1) = s.varl(j1) + lamExtra(count)*(-tlb(i1)+L2S(i1)/tdb(i1));
            s.varl(j1) = s.varl(j1) + lamExtra(count+1)*(tub(i1)-L2S(i1))/tdb(i1);
         else
            j1 = j1 - n_variable;
            s.linu(i1) = s.linu(i1) + lamExtra(count+3)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(i1) = s.linu(i1) + lamExtra(count+2)*(-tlb(j1)+L2S(j1))/tdb(j1);
            s.linl(i1) = s.linl(i1) + lamExtra(count)*(-tlb(j1)+L2S(j1)/tdb(j1));
            s.linl(i1) = s.linl(i1) + lamExtra(count+1)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(j1) = s.linu(j1) + lamExtra(count+3)*(tub(i1)-L2S(i1))/tdb(i1);
            s.linu(j1) = s.linu(j1) + lamExtra(count+2)*(-tlb(i1)+L2S(i1))/tdb(i1);
            s.linl(j1) = s.linl(j1) + lamExtra(count)*(-tlb(i1)+L2S(i1)/tdb(i1));
            s.linl(j1) = s.linl(j1) + lamExtra(count+1)*(tub(i1)-L2S(i1))/tdb(i1);
         end
      end
      count = count+4;
   end
end







   function [A,b,c,K] = setSedumi()
      Alam = [-speye(n_opt-1), spalloc(n_opt-1,1,0)];
      clam = spalloc(n_opt-1,1,0);
      As = zeros((n_variable+1)^2,n_opt);
      model0 = QOIobj;
      model0Var = model0.VarNames;
      varName = obj.VarNames;
      [~,id1,id2] = intersect(varName,model0Var);
      id1 = [1;id1+1];
      id2 = [1;id2+1];
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         N0 = model0.NormalizedNumerator;
         D0 = model0.NormalizedDenominator;
         targetN = zeros(n_variable+1,n_variable+1);
         targetD = zeros(n_variable+1,n_variable+1);
         for j = 1:length(id2)
            for k = 1:length(id2)
               targetN(id1(j),id1(k)) = N0(id2(j),id2(k));
               targetD(id1(j),id1(k)) = D0(id2(j),id2(k));
            end
         end
         targetMatrix = targetN;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         Coef0 = model0.NormalizedCoefMatrix;
         targetMatrix = zeros(n_variable+1,n_variable+1);
         for j = 1:length(id2)
            for k = 1:length(id2)
               targetMatrix(id1(j),id1(k)) = Coef0(id2(j),id2(k));
            end
         end
      end
      cs = targetMatrix(:);
      for i = 1:n_units
         As(:,i) = -Qunits{2*i-1}(:);
         As(:,i+n_units) = -Qunits{2*i}(:);
      end
      for i = 1:n_variable
         As(:,2*n_units+i) = -Qx{2*i-1}(:);
         As(:,2*n_units+n_variable+i) = -Qx{2*i}(:);
      end
      if n_extra ~= 0
         count = 1;
         for k = 1:n_extra
            As(:,2*(n_units+n_variable)+count) = -Qextra{count}(:);
            As(:,2*(n_units+n_variable)+count+1) = -Qextra{count+1}(:);
            As(:,2*(n_units+n_variable)+count+2) = -Qextra{count+2}(:);
            As(:,2*(n_units+n_variable)+count+3) = -Qextra{count+3}(:);
            count = count+4;
         end
      end          
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         As(:,end) = targetD;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         As(1,end) = 1;
      end
      At = [Alam; As];
      c = [clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = 1;
      K.l = n_opt-1;
      K.s = n_variable+1;
      A = At';    
   end

=======
function [minout,s] = sedumiminouterbound(obj, QOIobj, frac,abE)
% Calculate the outer bound of minimum values of the dataset
% target QOI subject to the constraints within the dataset through sedumi
% QOIobj - A B2BDC.B2Bmodel.Model object specifies the QOI
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015    Wenyu Li
%  Modified: Feb 16, 2016  Wenyu Li (Normalized sensitivity calculated)

bds = obj.calBound;
bds(:,1) = bds(:,1) - abE;
bds(:,2) = bds(:,2) + abE;
d = bds(:,2)-bds(:,1);
n_units = obj.DatasetUnits.Length;
units = obj.DatasetUnits.Values;
for n1 = 1:n_units
   my = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.dy;
   bds(n1,:) = (bds(n1,:) - my) / dy;
   d(n1) = d(n1)/dy;
end
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
n_variable = obj.Variables.Length;
vd = 2*ones(n_variable,1);
[Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
s = [];
nL = 0.5*length(Qx) - n_variable;
n_opt = 1+2*(n_units+n_variable+nL+2*n_extra);
if nL > 0
   tlb = obj.Variables.ExtraLinConstraint.LB;
   tub = obj.Variables.ExtraLinConstraint.UB;
   tdb = tub-tlb;
   L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
end

[A,b,c,K] = setSedumi;
pars.fid = 0;
[xopt,yopt,info] = sedumi(A,b,c,K,pars);
S = reshape(xopt(n_opt:end),n_variable+1,n_variable+1);
S = 0.5*(S+S');

if info.pinf ~= 0 || info.dinf ~= 0
   disp('Not both primal/dual problem are feasible, please review the results with cautions')
end
minout = yopt(end);
lamEU = -yopt(1:n_units);
lamEL = -yopt(n_units+1:2*n_units);
lamVU = -yopt(2*n_units+1:2*n_units+n_variable);
lamVL = -yopt(2*n_units+n_variable+1:2*(n_units+n_variable));
if nL > 0
   lamLU = -yopt(2*(n_units+n_variable)+1:2*(n_units+n_variable)+nL);
   lamLL = -yopt(2*(n_units+n_variable)+nL+1:2*(n_units+n_variable+nL));
end
if n_extra ~= 0
   lamExtra = -yopt(2*(n_units+nL+n_variable)+1:end-1);
end
s.expu = zeros(n_units,1);
s.expl = zeros(n_units,1);
s.expu(~idRQ) = lamEU(~idRQ).*d(~idRQ)*S(1,1);
s.expl(~idRQ) = lamEL(~idRQ).*d(~idRQ)*S(1,1);
for j1 = 1:sum(idRQ)
   ida = find(idRQ,j1);
   idU = ida(end);
   D = units(idU).SurrogateModel.NormalizedDenominator;
   s.expu(idU) = lamEU(idU)*d(idU)*dot(S,D);
   s.expl(idU) = lamEL(idU)*d(idU)*dot(S,D);
end
s.varu = (1.1*S(1,1) + S(2:end,1)).*lamVU.*vd;
s.varl = (1.1*S(1,1) - S(2:end,1)).*lamVL.*vd;
if nL > 0
   L2S = L2*S(2:end,1);
   s.linu = lamLU.*((-tlb-0.05*tdb)*S(1,1) + L2S).*tdb;
   s.linl = lamLL.*((tub+0.05*tdb)*S(1,1) - L2S).*tdb;
else
   s.linu = [];
   s.linl = [];
end
if n_extra ~= 0
   count = 1;
   for k1 = 1:n_extra
      i1 = extraIdx(k1,1);
      j1 = extraIdx(k1,2);
      if i1 <= n_variable
         if j1 <= n_variable
            s.varu(i1) = s.varu(i1) + 0.5*lamExtra(count+3)*(1-S(1,j1));
            s.varu(i1) = s.varu(i1) + 0.5*lamExtra(count+2)*(1+S(1,j1));
            s.varl(i1) = s.varl(i1) + 0.5*lamExtra(count)*(1+S(1,j1));
            s.varl(i1) = s.varl(i1) + 0.5*lamExtra(count+1)*(1-S(1,j1));
            s.varu(j1) = s.varu(j1) + 0.5*lamExtra(count+3)*(1-S(1,i1));
            s.varu(j1) = s.varu(j1) + 0.5*lamExtra(count+2)*(1+S(1,i1));
            s.varl(j1) = s.varl(j1) + 0.5*lamExtra(count)*(1+S(1,i1));
            s.varl(j1) = s.varl(j1) + 0.5*lamExtra(count+1)*(1-S(1,i1));
         else
            j1 = j1 - n_variable;
            s.varu(i1) = s.varu(i1) + lamExtra(count+3)*(tub(j1)-L2S(j1))/tdb(j1);
            s.varu(i1) = s.varu(i1) + lamExtra(count+2)*(-tlb(j1)+L2S(j1))/tdb(j1);
            s.varl(i1) = s.varl(i1) + lamExtra(count)*(-tlb(j1)+L2S(j1)/tdb(j1));
            s.varl(i1) = s.varl(i1) + lamExtra(count+1)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(j1) = s.linu(j1) + 0.5*lamExtra(count+3)*(1-S(1,i1));
            s.linu(j1) = s.linu(j1) + 0.5*lamExtra(count+2)*(1+S(1,i1));
            s.linl(j1) = s.linl(j1) + 0.5*lamExtra(count)*(1+S(1,i1));
            s.linl(j1) = s.linl(j1) + 0.5*lamExtra(count+1)*(1-S(1,i1));
         end
      else
         i1 = i1 - n_variable;
         if j1 <= n_variable
            s.linu(i1) = s.linu(i1) + 0.5*lamExtra(count+3)*(1-S(1,j1));
            s.linu(i1) = s.linu(i1) + 0.5*lamExtra(count+2)*(1+S(1,j1));
            s.linl(i1) = s.linl(i1) + 0.5*lamExtra(count)*(1+S(1,j1));
            s.linl(i1) = s.linl(i1) + 0.5*lamExtra(count+1)*(1-S(1,j1));
            s.varu(j1) = s.varu(j1) + lamExtra(count+3)*(tub(i1)-L2S(i1))/tdb(i1);
            s.varu(j1) = s.varu(j1) + lamExtra(count+2)*(-tlb(i1)+L2S(i1))/tdb(i1);
            s.varl(j1) = s.varl(j1) + lamExtra(count)*(-tlb(i1)+L2S(i1)/tdb(i1));
            s.varl(j1) = s.varl(j1) + lamExtra(count+1)*(tub(i1)-L2S(i1))/tdb(i1);
         else
            j1 = j1 - n_variable;
            s.linu(i1) = s.linu(i1) + lamExtra(count+3)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(i1) = s.linu(i1) + lamExtra(count+2)*(-tlb(j1)+L2S(j1))/tdb(j1);
            s.linl(i1) = s.linl(i1) + lamExtra(count)*(-tlb(j1)+L2S(j1)/tdb(j1));
            s.linl(i1) = s.linl(i1) + lamExtra(count+1)*(tub(j1)-L2S(j1))/tdb(j1);
            s.linu(j1) = s.linu(j1) + lamExtra(count+3)*(tub(i1)-L2S(i1))/tdb(i1);
            s.linu(j1) = s.linu(j1) + lamExtra(count+2)*(-tlb(i1)+L2S(i1))/tdb(i1);
            s.linl(j1) = s.linl(j1) + lamExtra(count)*(-tlb(i1)+L2S(i1)/tdb(i1));
            s.linl(j1) = s.linl(j1) + lamExtra(count+1)*(tub(i1)-L2S(i1))/tdb(i1);
         end
      end
      count = count+4;
   end
end







   function [A,b,c,K] = setSedumi()
      Alam = [-speye(n_opt-1), spalloc(n_opt-1,1,0)];
      clam = spalloc(n_opt-1,1,0);
      As = zeros((n_variable+1)^2,n_opt);
      model0 = QOIobj;
      model0Var = model0.VarNames;
      varName = obj.VarNames;
      [~,id1,id2] = intersect(varName,model0Var);
      id1 = [1;id1+1];
      id2 = [1;id2+1];
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         N0 = model0.NormalizedNumerator;
         D0 = model0.NormalizedDenominator;
         targetN = zeros(n_variable+1,n_variable+1);
         targetD = zeros(n_variable+1,n_variable+1);
         for j = 1:length(id2)
            for k = 1:length(id2)
               targetN(id1(j),id1(k)) = N0(id2(j),id2(k));
               targetD(id1(j),id1(k)) = D0(id2(j),id2(k));
            end
         end
         targetMatrix = targetN;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         Coef0 = model0.NormalizedCoefMatrix;
         targetMatrix = zeros(n_variable+1,n_variable+1);
         for j = 1:length(id2)
            for k = 1:length(id2)
               targetMatrix(id1(j),id1(k)) = Coef0(id2(j),id2(k));
            end
         end
      end
      cs = targetMatrix(:);
      for i = 1:n_units
         As(:,i) = -Qunits{2*i-1}(:);
         As(:,i+n_units) = -Qunits{2*i}(:);
      end
      for i = 1:n_variable
         As(:,2*n_units+i) = -Qx{2*i-1}(:);
         As(:,2*n_units+n_variable+i) = -Qx{2*i}(:);
      end
      if n_extra ~= 0
         count = 1;
         for k = 1:n_extra
            As(:,2*(n_units+n_variable)+count) = -Qextra{count}(:);
            As(:,2*(n_units+n_variable)+count+1) = -Qextra{count+1}(:);
            As(:,2*(n_units+n_variable)+count+2) = -Qextra{count+2}(:);
            As(:,2*(n_units+n_variable)+count+3) = -Qextra{count+3}(:);
            count = count+4;
         end
      end          
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         As(:,end) = targetD;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         As(1,end) = 1;
      end
      At = [Alam; As];
      c = [clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = 1;
      K.l = n_opt-1;
      K.s = n_variable+1;
      A = At';    
   end

>>>>>>> Stashed changes
end