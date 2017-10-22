<<<<<<< Updated upstream
function [y,s] = sedumiconsisrel(obj,yin,opt,abE)
% Calculate outer bound for the consistency measure of a
% B2BDC.B2Bdataset.Dateset object through sedumi
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015    Wenyu Li
%  Modified: September 1, 2015  Wenyu Li (Normalized sensitivity calculated) 

frac = opt.ExtraLinFraction;
n_units = obj.DatasetUnits.Length;
n_variable = obj.Variables.Length;
units = obj.DatasetUnits.Values;
l = [obj.DatasetUnits.Values.LowerBound]';
u = [obj.DatasetUnits.Values.UpperBound]';
l = l-abE;
u = u+abE;
d = [obj.DatasetUnits.Values.ObservedValue]';
for j = 1:n_units
   my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
   l(j) = (l(j)-my)/dy;
   u(j) = (u(j)-my)/dy;
   d(j) = (d(j)-my)/dy;
end
qbd = u-l;
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
% vd = vu-vl;
vd = 2*ones(n_variable,1);
t_1 = yin;
t_2 = 1;
s = [];
while t_2-t_1 > opt.TolConsis
   t_m = 0.5*(t_1+t_2);
   bds = [(1-t_m)*l+t_m*d , (1-t_m)*u+t_m*d];
   [Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
   nL = 0.5*length(Qx) - n_variable;
   if nL > 0
      tlb = obj.Variables.ExtraLinConstraint.LB;
      tub = obj.Variables.ExtraLinConstraint.UB;
      tdb = tub-tlb;
      L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
   end
   n_opt = 2*(n_units+n_variable+nL+2*n_extra)+1;
   [A,b,c,K] = setSedumi;
   pars.fid = 0;
   pars.eps = 1e-7;
   [xopt,yopt] = sedumi(A,b,c,K,pars);
   if yopt(end) >= 0
      S = reshape(xopt(n_opt+1:end),n_variable+1,n_variable+1);
      S = 0.5*(S+S');
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
      t_2 = t_m;
      s.expu = zeros(n_units,1);
      s.expl = zeros(n_units,1);
      s.expu(~idRQ) = lamEU(~idRQ).*qbd(~idRQ)*S(1,1);
      s.expl(~idRQ) = lamEL(~idRQ).*qbd(~idRQ)*S(1,1);
      for j = 1:sum(idRQ)
         id1 = find(idRQ,j);
         idU = id1(end);
         D = units(idU).SurrogateModel.NormalizedDenominator;
         s.expu(idU) = lamEU(idU)*qbd(idU)*dot(S,D);
         s.expl(idU) = lamEL(idU)*qbd(idU)*dot(S,D);
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
      s.expu(abs(s.expu)<1e-5) = 0;
      s.expl(abs(s.expl)<1e-5) = 0;
      s.varu(abs(s.varu)<1e-5) = 0;
      s.varl(abs(s.varl)<1e-5) = 0;
      s.linu(abs(s.linu)<1e-5) = 0;
      s.linl(abs(s.linl)<1e-5) = 0;
   else
      t_1 = t_m;
   end
end
y = 0.5*(t_2+t_1);
   
   
   
   
   
   
   
   

 




   function [A,b,c,K] = setSedumi()
      Ac = [-ones(1,n_opt-1), 0];
      cc = -1;
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
      As(:,end) = vec(speye(n_variable+1));
      cs = spalloc((n_variable+1)^2,1,0);
      At = [Ac; Alam; As];
      c = [cc; clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = 1;
      K.l = n_opt;
      K.s = n_variable+1;
      A = At';
   end

=======
function [y,s] = sedumiconsisrel(obj,yin,opt,abE)
% Calculate outer bound for the consistency measure of a
% B2BDC.B2Bdataset.Dateset object through sedumi
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015    Wenyu Li
%  Modified: September 1, 2015  Wenyu Li (Normalized sensitivity calculated) 

frac = opt.ExtraLinFraction;
n_units = obj.DatasetUnits.Length;
n_variable = obj.Variables.Length;
units = obj.DatasetUnits.Values;
l = [obj.DatasetUnits.Values.LowerBound]';
u = [obj.DatasetUnits.Values.UpperBound]';
l = l-abE;
u = u+abE;
d = [obj.DatasetUnits.Values.ObservedValue]';
for j = 1:n_units
   my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
   l(j) = (l(j)-my)/dy;
   u(j) = (u(j)-my)/dy;
   d(j) = (d(j)-my)/dy;
end
qbd = u-l;
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
% vd = vu-vl;
vd = 2*ones(n_variable,1);
t_1 = yin;
t_2 = 1;
s = [];
while t_2-t_1 > opt.TolConsis
   t_m = 0.5*(t_1+t_2);
   bds = [(1-t_m)*l+t_m*d , (1-t_m)*u+t_m*d];
   [Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
   nL = 0.5*length(Qx) - n_variable;
   if nL > 0
      tlb = obj.Variables.ExtraLinConstraint.LB;
      tub = obj.Variables.ExtraLinConstraint.UB;
      tdb = tub-tlb;
      L2 = L(n_variable+1:end,:) .* repmat(tdb,1,n_variable);
   end
   n_opt = 2*(n_units+n_variable+nL+2*n_extra)+1;
   [A,b,c,K] = setSedumi;
   pars.fid = 0;
   pars.eps = 1e-7;
   [xopt,yopt] = sedumi(A,b,c,K,pars);
   if yopt(end) >= 0
      S = reshape(xopt(n_opt+1:end),n_variable+1,n_variable+1);
      S = 0.5*(S+S');
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
      t_2 = t_m;
      s.expu = zeros(n_units,1);
      s.expl = zeros(n_units,1);
      s.expu(~idRQ) = lamEU(~idRQ).*qbd(~idRQ)*S(1,1);
      s.expl(~idRQ) = lamEL(~idRQ).*qbd(~idRQ)*S(1,1);
      for j = 1:sum(idRQ)
         id1 = find(idRQ,j);
         idU = id1(end);
         D = units(idU).SurrogateModel.NormalizedDenominator;
         s.expu(idU) = lamEU(idU)*qbd(idU)*dot(S,D);
         s.expl(idU) = lamEL(idU)*qbd(idU)*dot(S,D);
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
      s.expu(abs(s.expu)<1e-5) = 0;
      s.expl(abs(s.expl)<1e-5) = 0;
      s.varu(abs(s.varu)<1e-5) = 0;
      s.varl(abs(s.varl)<1e-5) = 0;
      s.linu(abs(s.linu)<1e-5) = 0;
      s.linl(abs(s.linl)<1e-5) = 0;
   else
      t_1 = t_m;
   end
end
y = 0.5*(t_2+t_1);
   
   
   
   
   
   
   
   

 




   function [A,b,c,K] = setSedumi()
      Ac = [-ones(1,n_opt-1), 0];
      cc = -1;
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
      As(:,end) = vec(speye(n_variable+1));
      cs = spalloc((n_variable+1)^2,1,0);
      At = [Ac; Alam; As];
      c = [cc; clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = 1;
      K.l = n_opt;
      K.s = n_variable+1;
      A = At';
   end

>>>>>>> Stashed changes
end