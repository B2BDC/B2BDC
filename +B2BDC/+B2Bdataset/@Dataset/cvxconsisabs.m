function [y,s] = cvxconsisabs(obj,yin,frac,tolerance)
% Calculate outer bound for the consistency measure of a
% B2BDC.B2Bdataset.Dateset object through cvx
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included

% Created: June 15, 2015     Wenyu Li
%  Modified: June 25, 2015     Wenyu Li
%  Modified: July 2, 2015     Wenyu Li  (sensitivity added)
%  Modified July 12, 2015    Wenyu Li (higher order calculation added)

l = [obj.DatasetUnits.Values.LowerBound]';
u = [obj.DatasetUnits.Values.UpperBound]';
dmin = min(u-l);
n_units = obj.DatasetUnits.Length;
n_variable = obj.Variables.Length;
t_1 = yin;
t_2 = 0.51*(dmin);
s = [];
while t_2-t_1 > tolerance
   t_m = 0.5*(t_1+t_2);
   bds = [l+t_m , u-t_m];
   [Qunits, Qx, Qextra, n_extra]  = obj.getInequalQuad(bds,frac);
   cvx_begin sdp quiet
   cvx_precision('best');
   variable lamlb(n_units)
   variable lamub(n_units)
   variable lamxl(n_variable)
   variable lamxu(n_variable)
   if n_extra ~= 0
      variable lamextra(4*n_extra)
   end
   cvxconstr = zeros(n_variable+1,n_variable+1);
   for i = 1:n_units
      cvxconstr = cvxconstr + lamub(i)*Qunits{2*i-1}+...
         lamlb(i)*Qunits{2*i};
   end
   for i = 1:n_variable
      cvxconstr = cvxconstr + lamxu(i)*Qx{2*i-1}+...
         lamxl(i)*Qx{2*i};
   end
   if n_extra ~= 0
      count = 1;
      for k = 1:n_extra
         cvxconstr = cvxconstr + lamextra(count)*Qextra{count}+...
            lamextra(count+1)*Qextra{count+1} + lamextra(count+2)*Qextra{count+2}+...
            lamextra(count+3)*Qextra{count+3};
         count = count+4;
      end
   end
   if n_extra ~= 0
      minimize sum(lamlb)+sum(lamub)+sum(lamxl)+sum(lamxu)+sum(lamextra)
   else
      minimize sum(lamlb)+sum(lamub)+sum(lamxl)+sum(lamxu)
   end
   subject to
   cvxconstr >= zeros(n_variable+1,n_variable+1);
   lamlb >= 0;
   lamub >= 0;
   lamxl >= 0;
   lamxu >= 0;
   if n_extra ~= 0
      lamextra >= 0;
   end
   if n_extra ~= 0
       sum(lamlb)+sum(lamub)+sum(lamxl)+sum(lamxu)+sum(lamextra) >= 1;
   else
       sum(lamlb)+sum(lamub)+sum(lamxl)+sum(lamxu) >= 1;
   end
   cvx_end
   if strcmp(cvx_status,'Solved')
      t_2 = t_m;
      s.expl = -lamlb;
      s.expu = lamub;
      s.varl = -lamxl;
      s.varu = lamxu;
   else
      t_1 = t_m;
   end
end
y = 0.5*(t_2+t_1);