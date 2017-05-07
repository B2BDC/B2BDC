function [y,s] = sedumiconsisquadrel(obj,opt,abE)
% Calculate outer bound for consistency measure (in relative value) of a
% B2BDC.B2Bdataset.Dataset object through cvx for the case when all
% surrogate models are quadratic functions.
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included
% The problem is formed in the dual format, and the dual variables are
% arranged as following:
%           y = [ lamub, lamlb, lamxu, lamxl, lamextra, ro ]
%  lamub  -  Lagrangian multiplier wrt dataset units upper bound
%  lamul  -  Lagrangian multiplier wrt dataset units lower bound
%  lamxu  -  Lagrangian multiplier wrt variable upper bound
%  lamxl  -  Lagrangian multiplier wrt variable lower bound
% lamextra - Lagrangian multiplier wrt extra linear combination pairs
%    ro   -  optimizing variable

% Created July 30, 2015    Wenyu Li
%  Modified: Feb 16, 2016  Wenyu Li (Normalized sensitivity calculated) 

n_units = obj.DatasetUnits.Length;
frac = opt.ExtraLinFraction;
unitsUB = [obj.DatasetUnits.Values.UpperBound]' - [obj.DatasetUnits.Values.ObservedValue]';
unitsLB = [obj.DatasetUnits.Values.LowerBound]' - [obj.DatasetUnits.Values.ObservedValue]';
unitsUB = unitsUB+abE;
unitsLB = unitsLB-abE;
n_variable = obj.Variables.Length;
bds = obj.calBound;
bds(:,1) = bds(:,1) - abE;
bds(:,2) = bds(:,2) + abE;
d = bds(:,2)-bds(:,1);
% for j = 1:n_units
%    my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
%    dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
%    bds(j,:) = (bds(j,:) - my) / dy;
%    d(j) = d(j)/dy;
%    unitsUB(j) = unitsUB(j) / dy;
%    unitsLB(j) = unitsLB(j) / dy;
% end
[Qunits, Qx, Qextra, n_extra, extraIdx, L]  = obj.getInequalQuad(bds,frac);
nL = length(Qx) - n_variable;
n_opt = 1+2*n_units+n_variable+nL+n_extra;
vu = [obj.Variables.Values.UpperBound]';
vl = [obj.Variables.Values.LowerBound]';
% vu = 0.5*ones(n_variable,1);
% vl = -0.5*vu;
% vd = 2*ones(n_variable,1);
vd = vu - vl;
if nL > 0
   tlb = obj.Variables.ExtraLinConstraint.LB;
   tub = obj.Variables.ExtraLinConstraint.UB;
   tdb = tub-tlb;
   L2 = L(n_variable+1:end,:);
end

[A,b,c,K] = setSedumi;
pars.fid = 0;
[xopt,yopt,info] = sedumi(A,b,c,K,pars);
S = reshape(xopt(n_opt+2:end),n_variable+1,n_variable+1);
S = 0.5*(S+S');
if info.pinf ~= 0 || info.dinf ~= 0 || info.numerr ~= 0
   disp('Not both primal/dual problem are feasible or there exsit')
   disp('numerical inaccuracy, please review the results with cautions')
end
y = yopt(end);
lamEU = yopt(1:n_units);
lamEL = yopt(n_units+1:2*n_units);
lamV = yopt(2*n_units+1:2*n_units+n_variable);
if nL > 0
   lamL = yopt(2*n_units+n_variable+1:2*n_units+n_variable+nL);
end
if n_extra ~= 0
   lamExtra = yopt(2*n_units+n_variable+nL+1:end-1);
end

s.expu = lamEU.*d*S(1,1)*(1+xopt(1)-xopt(2));
s.expl = lamEL.*d*S(1,1)*(1+xopt(1)-xopt(2));
% s.varu = 1.05*lamVU.*vd.^2;
% s.varl = 1.05*lamVL.*vd.^2;
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

s.expu(abs(s.expu)<1e-5) = 0;
s.expl(abs(s.expl)<1e-5) = 0;
s.varu(abs(s.varu)<1e-5) = 0;
s.varl(abs(s.varl)<1e-5) = 0;
s.linu(abs(s.linu)<1e-5) = 0;
s.linl(abs(s.linl)<1e-5) = 0;




   function [A,b,c,K] = setSedumi()
      Ac = [unitsUB', -unitsLB', spalloc(1,n_opt-2*n_units,0);
         -unitsUB', unitsLB', spalloc(1,n_opt-2*n_units,0)];
      cc = [1;-1];
      Alam = [-speye(n_opt-1), spalloc(n_opt-1,1,0)];
      clam = spalloc(n_opt-1,1,0);
      As = zeros((n_variable+1)^2,n_opt);
      for i = 1:n_units
         As(:,i) = -Qunits{2*i-1}(:);
         As(:,i+n_units) = -Qunits{2*i}(:);
      end
      for i = 1:n_variable
         As(:,2*n_units+i) = -Qx{i}(:);
      end
      for i = 1:nL
         As(:,2*n_units+n_variable+i) = -Qx{n_variable+i}(:);
      end
      if n_extra ~= 0
         for k = 1:n_extra
            As(:,2*n_units+n_variable+nL+k) = -Qextra{k}(:);
         end
      end           
      As(1,end) = -1;
      cs = spalloc((n_variable+1)^2,1,0);
      At = [Ac; Alam; As];
      c = [cc; clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = -1;
      K.l = 2*n_units+n_variable+nL+n_extra+2;
      K.s = n_variable+1;
      A = At';    
   end

end