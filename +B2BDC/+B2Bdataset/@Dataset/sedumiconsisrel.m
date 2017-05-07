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
name2 = obj.VarNames;
l = [obj.DatasetUnits.Values.LowerBound]';
u = [obj.DatasetUnits.Values.UpperBound]';
l = l-abE;
u = u+abE;
d = [obj.DatasetUnits.Values.ObservedValue]';
% for j = 1:n_units
%    my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
%    dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
%    l(j) = (l(j)-my)/dy;
%    u(j) = (u(j)-my)/dy;
%    d(j) = (d(j)-my)/dy;
% end
qbd = u-l;
vu = [obj.Variables.Values.UpperBound]';
vl = [obj.Variables.Values.LowerBound]';
vd = vu-vl;
% vd = 2*ones(n_variable,1);
t_1 = yin;
t_2 = 1;
s = [];
while t_2-t_1 > opt.TolConsis
   t_m = 0.5*(t_1+t_2);
   bds = [(1-t_m)*l+t_m*d , (1-t_m)*u+t_m*d];
   [Qunits, Qx, Qextra, n_extra, extraIdx,L,idRQ]  = obj.getInequalQuad(bds,frac);
   nL = length(Qx) - n_variable;
   if nL > 0
      tlb = obj.Variables.ExtraLinConstraint.LB;
      tub = obj.Variables.ExtraLinConstraint.UB;
      tdb = tub-tlb;
      L2 = L(n_variable+1:end,:);
   end
   n_opt = 2*n_units+n_variable+nL+n_extra+1;
   [A,b,c,K] = setSedumi;
   pars.fid = 0;
   pars.eps = 1e-7;
   [xopt,yopt] = sedumi(A,b,c,K,pars);
   if yopt(end) >= 0
%       S = reshape(xopt(n_opt+1:end),n_variable+1,n_variable+1);
%       S = 0.5*(S+S');
%       lamEU = yopt(1:n_units);
%       lamEL = yopt(n_units+1:2*n_units);
%       lamV = yopt(2*n_units+1:2*n_units+n_variable);
%       if nL > 0
%          lamL = yopt(2*n_units+n_variable+1:2*n_units+n_variable+nL);
%       end
%       if n_extra ~= 0
%          lamExtra = yopt(2*(n_units+n_variable+nL)+1:end-1);
%       end
      t_2 = t_m;
%       s.expu = zeros(n_units,1);
%       s.expl = zeros(n_units,1);
%       s.expu(~idRQ) = lamEU(~idRQ).*qbd(~idRQ)*S(1,1);
%       s.expl(~idRQ) = lamEL(~idRQ).*qbd(~idRQ)*S(1,1);
%       for j = 1:sum(idRQ)
%          id1 = find(idRQ,j);
%          idU = id1(end);
%          D = units(idU).SurrogateModel.Denominator;
%          name1 = units(idU).SurrogateModel.VarNames;
%          [~,~,id] = intersect(name1,name2,'stable');
%          id = [1;id+1];
%          s.expu(idU) = lamEU(idU)*qbd(idU)*trace(S(id,id)*D');
%          s.expl(idU) = lamEL(idU)*qbd(idU)*trace(S(id,id)*D');
%       end
%       s.varu = (-vl*S(1,1) + S(2:end,1)).*lamV.*vd;
%       s.varl = (vu*S(1,1) - S(2:end,1)).*lamV.*vd;
%       if nL > 0
%          L2S = L2*S(2:end,1);
%          s.linu = lamL.*(-tlb*S(1,1) + L2S).*tdb;
%          s.linl = lamL.*(tub*S(1,1) - L2S).*tdb;
%       else
%          s.linu = [];
%          s.linl = [];
%       end
      
%       if n_extra ~= 0
%          LS = L*S(2:end,1);
%          for k1 = 1:n_extra
%             i1 = extraIdx(k1,1);
%             j1 = extraIdx(k1,1);
%             c1 = extraIdx(k1,1);
%             if c1 == 1
%                if i1 <= n_variable
%                   li = vl(i1);
%                   if j1 <= n_variable
%                      lj = vl(j1);
%                      s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   else
%                      lj = tlb(j1-n_variable);
%                      s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   end
%                   s.varl(i1) = s.varl(i1) + vd(i1)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
%                else
%                   li = tlb(i1-n_variable);
%                   if j1 <= n_variable
%                      lj = vl(j1);
%                      s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   else
%                      lj = tlb(j1-n_variable);
%                      s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   end
%                   s.linl(i1-n_variable) = s.linl(i1-n_variable) + ...
%                      tdb(i1-n_variable)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
%                end
%             elseif c1 == 2
%                if i1 <= n_variable
%                   li = vl(i1);
%                   if j1 <= n_variable
%                      uj = vu(j1);
%                      s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   else
%                      uj = tub(j1-n_variable);
%                      s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   end
%                   s.varl(i1) = s.varl(i1) + vd(i1)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
%                else
%                   li = tlb(i1-n_variable);
%                   if j1 <= n_variable
%                      uj = vu(j1);
%                      s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   else
%                      uj = tub(j1-n_variable);
%                      s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(-li*S(1,1)+LS(i1));
%                   end
%                   s.linl(i1-n_variable) = s.linl(i1-n_variable) + ...
%                      tdb(i1-n_variable)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
%                end
%             elseif c1 == 3
%                if i1 <= n_variable
%                   ui = vu(i1);
%                   if j1 <= n_variable
%                      lj = vl(j1);
%                      s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   else
%                      lj = tlb(j1-n_variable);
%                      s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   end
%                   s.varu(i1) = s.varu(i1) + vd(i1)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
%                else
%                   ui = tub(i1-n_variable);
%                   if j1 <= n_variable
%                      lj = vl(j1);
%                      s.varl(j1) = s.varl(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   else
%                      lj = tlb(j1-n_variable);
%                      s.linl(j1-n_variable) = s.linl(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   end
%                   s.linu(i1-n_variable) = s.linu(i1-n_variable) + ...
%                      tdb(i1-n_variable)*lamExtra(k1)*(-lj*S(1,1)+LS(j1));
%                end
%             else
%                if i1 <= n_variable
%                   ui = vu(i1);
%                   if j1 <= n_variable
%                      uj = vu(j1);
%                      s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   else
%                      uj = tub(j1-n_variable);
%                      s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   end
%                   s.varu(i1) = s.varu(i1) + vd(i1)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
%                else
%                   ui = tub(i1-n_variable);
%                   if j1 <= n_variable
%                      uj = vu(j1);
%                      s.varu(j1) = s.varu(j1) + vd(j1)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   else
%                      uj = tub(j1-n_variable);
%                      s.linu(j1-n_variable) = s.linu(j1-n_variable) + ...
%                         tdb(j1-n_variable)*lamExtra(k1)*(ui*S(1,1)-LS(i1));
%                   end
%                   s.linu(i1-n_variable) = s.linu(i1-n_variable) + ...
%                      tdb(i1-n_variable)*lamExtra(k1)*(uj*S(1,1)-LS(j1));
%                end
%             end
%          end
%       end

      
%       s.expu(abs(s.expu)<1e-5) = 0;
%       s.expl(abs(s.expl)<1e-5) = 0;
%       s.varu(abs(s.varu)<1e-5) = 0;
%       s.varl(abs(s.varl)<1e-5) = 0;
%       s.linu(abs(s.linu)<1e-5) = 0;
%       s.linl(abs(s.linl)<1e-5) = 0;
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

end