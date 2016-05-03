function [y,s] = sedumiconsisabs(obj,opt,abE)
% Calculate outer bound for the consistency measure of a
% B2BDC.B2Bdataset.Dateset object through sedumi
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: August 7, 2015    Wenyu Li
%  Modified: September 1, 2015  Wenyu Li (Normalized sensitivity calculated) 

frac = opt.ExtraLinFraction;
l = [obj.DatasetUnits.Values.LowerBound]';
u = [obj.DatasetUnits.Values.UpperBound]';
l = l-abE;
u = u+abE;
d = u-l;
for j = 1:n_units
   my = obj.DatasetUnits.Values(j).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(j).SurrogateModel.yScale.dy;
   l(j) = (l(j)-my)/dy;
   u(j) = (u(j)-my)/dy;
   d(j) = d(j)/dy;
end
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
dmin = min(d);
n_units = obj.DatasetUnits.Length;
n_variable = obj.Variables.Length;
vu = ones(n_variable,1);
vl = -vu;
vd = 2*vu;
t_1 = yin;
t_2 = 0.51*(dmin);
s = [];
while t_2-t_1 > opt.TolConsis
   t_m = 0.5*(t_1+t_2);
   bds = [l+t_m , u-t_m];
   [Qunits, Qx, Qextra, n_extra, extraIdx]  = obj.getInequalQuad(bds,frac);
   n_opt = 2*(n_units+n_variable+2*n_extra)+1;
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
      if n_extra ~= 0
         lamExtra = yopt(2*(n_units+n_variable)+1:end-1);
      end
      t_2 = t_m;
      s.expu = lamEU.*d*S(1,1);
      s.expl = lamEL.*d*S(1,1);
      s.varu = 1.05*lamVU.*vd.^2;
      s.varl = 1.05*lamVL.*vd.^2;
      if n_extra ~= 0
         count = 1;
         for k1 = 1:n_extra
            i1 = extraIdx(k1,1);
            j1 = extraIdx(k1,2);
            tmpLam = lamExtra(count:count+3);
            s.varu(i1) = s.varu(i1)-vd(i1)*(tmpLam(3)*(vl(j1)*S(1,1)-S(1,j1))+tmpLam(4)*(-vu(j1)*S(1,1)+S(1,j1)));
            s.varu(j1) = s.varu(j1)-vd(j1)*(tmpLam(2)*(vl(i1)*S(1,1)-S(1,i1))+tmpLam(4)*(-vu(i1)*S(1,1)+S(1,i1)));
            s.varl(i1) = s.varl(i1)-vd(i1)*(tmpLam(1)*(vl(j1)*S(1,1)-S(1,j1))+tmpLam(2)*(-vu(j1)*S(1,1)+S(1,j1)));
            s.varl(j1) = s.varl(j1)-vd(j1)*(tmpLam(1)*(vl(i1)*S(1,1)-S(1,i1))+tmpLam(3)*(-vu(i1)*S(1,1)+S(1,i1)));
            count = count+4;
         end
      end
      s.expu(abs(s.expu)<1e-5) = 0;
      s.expl(abs(s.expl)<1e-5) = 0;
      s.varu(abs(s.varu)<1e-5) = 0;
      s.varl(abs(s.varl)<1e-5) = 0;
   else
      t_1 = t_m;
   end
end
y = 0.5*(t_2+t_1);
   
   
   
   
   
   
   
   






   function [A,b,c,K] = setSedumi()
      Ac = [-ones(1,n_opt-1),0];
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