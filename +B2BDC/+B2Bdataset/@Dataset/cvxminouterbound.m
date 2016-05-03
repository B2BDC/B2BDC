function [minout,minSensitivity] = cvxminouterbound(obj,QOIobj, frac)
% Calculate the outer bound of minimum values of the dataset
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

bds = [[obj.DatasetUnits.Values.LowerBound]' , [obj.DatasetUnits.Values.UpperBound]'];
n_units = obj.DatasetUnits.Length;
n_variable = obj.Variables.Length;
[Qunits, Qx, Qextra, n_extra]  = obj.getInequalQuad(bds,frac);
cvx_begin sdp quiet
variable tmin
variable lamlb(n_units)
variable lamub(n_units)
variable lamxl(n_variable)
variable lamxu(n_variable)
if n_extra ~= 0
   variable lamextra(4*n_extra)
end
model0 = QOIobj;
model0Var = model0.VarNames;
varName = obj.VarNames;
[~,id1,id2] = intersect(varName,model0Var);
id1 = [1;id1+1];
id2 = [1;id2+1];
if isa(model0,'B2BDC.B2Bmodels.RQModel')
   N0 = model0.Numerator;
   D0 = model0.Denominator;
   targetN = zeros(n_variable+1,n_variable+1);
   targetD = zeros(n_variable+1,n_variable+1);
   for j = 1:length(id2)
      for k = 1:length(id2)
         targetN(id1(j),id1(k)) = N0(id2(j),id2(k));
         targetD(id1(j),id1(k)) = D0(id2(j),id2(k));
      end
   end
   cvxconstr = targetN-tmin*targetD;
elseif isa(model0,'B2BDC.B2Bmodels.QModel')
   Coef0 = model0.CoefMatrix;
   targetMatrix1 = zeros(n_variable+1,n_variable+1);
   targetMatrix2 = zeros(n_variable+1,n_variable+1);
   targetMatrix2(1,1) = 1;
   for j = 1:length(id2)
      for k = 1:length(id2)
         targetMatrix1(id1(j),id1(k)) = Coef0(id2(j),id2(k));
      end
   end
   cvxconstr = targetMatrix1-tmin*targetMatrix2;
end
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
maximize tmin
subject to
cvxconstr >= zeros(n_variable+1,n_variable+1);
lamlb >= 0;
lamub >= 0;
lamxl >= 0;
lamxu >= 0;
if n_extra ~= 0
   lamextra >= 0;
end
cvx_end
minout = cvx_optval;
minSensitivity.expl = lamlb;
minSensitivity.expu = -lamub;
minSensitivity.varl = lamxl;
minSensitivity.varu = -lamxu;
