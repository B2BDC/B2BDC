<<<<<<< Updated upstream
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
for j = 1:n_units
   my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
   bds(j,:) = (bds(j,:) - my) / dy;
   d(j) = d(j)/dy;
   unitsUB(j) = unitsUB(j) / dy;
   unitsLB(j) = unitsLB(j) / dy;
end
[Qunits, Qx, Qextra, n_extra, extraIdx, L]  = obj.getInequalQuad_old(bds,frac);
nL = 0.5*length(Qx) - n_variable;
n_opt = 1+2*(n_units+n_variable+nL+2*n_extra);
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
% vu = 0.5*ones(n_variable,1);
% vl = -0.5*vu;
vd = 2*ones(n_variable,1);
if nL > 0
   tlb = obj.Variables.ExtraLinConstraint.LB;
   tub = obj.Variables.ExtraLinConstraint.UB;
   tdb = tub-tlb;
%    tlb = tlb./tdb;
%    tub = tub./tdb;
%    vu = [tub;vu];
%    vl = [tlb;vl];
   L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
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
lamVU = yopt(2*n_units+1:2*n_units+n_variable);
lamVL = yopt(2*n_units+n_variable+1:2*(n_units+n_variable));
if nL > 0
   lamLU = yopt(2*(n_units+n_variable)+1:2*(n_units+n_variable)+nL);
   lamLL = yopt(2*(n_units+n_variable)+nL+1:2*(n_units+n_variable+nL));
end
if n_extra ~= 0
   lamExtra = yopt(2*(n_units+n_variable+nL)+1:end-1);
end

s.expu = lamEU.*d*S(1,1)*(1+xopt(1)-xopt(2));
s.expl = lamEL.*d*S(1,1)*(1+xopt(1)-xopt(2));
% s.varu = 1.05*lamVU.*vd.^2;
% s.varl = 1.05*lamVL.*vd.^2;
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
         As(:,2*n_units+i) = -Qx{2*i-1}(:);
         As(:,2*n_units+n_variable+i) = -Qx{2*i}(:);
      end
      for i = 1:nL
         As(:,2*(n_units+n_variable)+i) = -Qx{2*(n_variable+i)-1}(:);
         As(:,2*(n_units+n_variable)+nL+i) = -Qx{2*(n_variable+i)}(:);
      end
      if n_extra ~= 0
         count = 1;
         for k = 1:n_extra
            As(:,2*(n_units+n_variable+nL)+count) = -Qextra{count}(:);
            As(:,2*(n_units+n_variable+nL)+count+1) = -Qextra{count+1}(:);
            As(:,2*(n_units+n_variable+nL)+count+2) = -Qextra{count+2}(:);
            As(:,2*(n_units+n_variable+nL)+count+3) = -Qextra{count+3}(:);
            count = count+4;
         end
      end           
      As(1,end) = -1;
      cs = spalloc((n_variable+1)^2,1,0);
      At = [Ac; Alam; As];
      c = [cc; clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = -1;
      K.l = 2*(n_units+n_variable+nL+2*n_extra+1);
      K.s = n_variable+1;
      A = At';    
   end

=======
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
for j = 1:n_units
   my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
   bds(j,:) = (bds(j,:) - my) / dy;
   d(j) = d(j)/dy;
   unitsUB(j) = unitsUB(j) / dy;
   unitsLB(j) = unitsLB(j) / dy;
end
[Qunits, Qx, Qextra, n_extra, extraIdx, L]  = obj.getInequalQuad_old(bds,frac);
nL = 0.5*length(Qx) - n_variable;
n_opt = 1+2*(n_units+n_variable+nL+2*n_extra);
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
% vu = 0.5*ones(n_variable,1);
% vl = -0.5*vu;
vd = 2*ones(n_variable,1);
if nL > 0
   tlb = obj.Variables.ExtraLinConstraint.LB;
   tub = obj.Variables.ExtraLinConstraint.UB;
   tdb = tub-tlb;
%    tlb = tlb./tdb;
%    tub = tub./tdb;
%    vu = [tub;vu];
%    vl = [tlb;vl];
   L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
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
lamVU = yopt(2*n_units+1:2*n_units+n_variable);
lamVL = yopt(2*n_units+n_variable+1:2*(n_units+n_variable));
if nL > 0
   lamLU = yopt(2*(n_units+n_variable)+1:2*(n_units+n_variable)+nL);
   lamLL = yopt(2*(n_units+n_variable)+nL+1:2*(n_units+n_variable+nL));
end
if n_extra ~= 0
   lamExtra = yopt(2*(n_units+n_variable+nL)+1:end-1);
end

s.expu = lamEU.*d*S(1,1)*(1+xopt(1)-xopt(2));
s.expl = lamEL.*d*S(1,1)*(1+xopt(1)-xopt(2));
% s.varu = 1.05*lamVU.*vd.^2;
% s.varl = 1.05*lamVL.*vd.^2;
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
         As(:,2*n_units+i) = -Qx{2*i-1}(:);
         As(:,2*n_units+n_variable+i) = -Qx{2*i}(:);
      end
      for i = 1:nL
         As(:,2*(n_units+n_variable)+i) = -Qx{2*(n_variable+i)-1}(:);
         As(:,2*(n_units+n_variable)+nL+i) = -Qx{2*(n_variable+i)}(:);
      end
      if n_extra ~= 0
         count = 1;
         for k = 1:n_extra
            As(:,2*(n_units+n_variable+nL)+count) = -Qextra{count}(:);
            As(:,2*(n_units+n_variable+nL)+count+1) = -Qextra{count+1}(:);
            As(:,2*(n_units+n_variable+nL)+count+2) = -Qextra{count+2}(:);
            As(:,2*(n_units+n_variable+nL)+count+3) = -Qextra{count+3}(:);
            count = count+4;
         end
      end           
      As(1,end) = -1;
      cs = spalloc((n_variable+1)^2,1,0);
      At = [Ac; Alam; As];
      c = [cc; clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = -1;
      K.l = 2*(n_units+n_variable+nL+2*n_extra+1);
      K.s = n_variable+1;
      A = At';    
   end

>>>>>>> Stashed changes
end