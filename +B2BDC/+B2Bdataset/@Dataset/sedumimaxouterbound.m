function [maxout,maxSensitivity] = sedumimaxouterbound(obj,QOIobj, frac, abE)
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
d = bds(:,2)-bds(:,1);
n_units = obj.DatasetUnits.Length;
for n1 = 1:n_units
   my = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.my;
   dy = obj.DatasetUnits.Values(n1).SurrogateModel.yScale.dy;
   bds(n1,:) = (bds(n1,:) - my) / dy;
   d(n1) = d(n1)/dy;
end
% vu = [obj.Variables.Values.UpperBound]';
% vl = [obj.Variables.Values.LowerBound]';
n_variable = obj.Variables.Length;
vu = ones(n_variable,1);
vl = -vu;
vd = 2*ones(n_variable,1);
[Qunits, Qx, Qextra, n_extra, extraIdx]  = obj.getInequalQuad(bds,frac);
maxSensitivity = [];
n_opt = 1+2*(n_units+n_variable+2*n_extra);

[A,b,c,K] = setSedumi;
pars.fid = 0;
[xopt,yopt,info] = sedumi(A,b,c,K,pars);
S = reshape(xopt(n_opt:end),n_variable+1,n_variable+1);
S = 0.5*(S+S');

if info.pinf ~= 0 || info.dinf ~= 0
   disp('Not both primal/dual problem are feasible, please review the results with cautions')
end
maxout = yopt(end);
lamEU = yopt(1:n_units);
lamEL = yopt(n_units+1:2*n_units);
lamVU = yopt(2*n_units+1:2*n_units+n_variable);
lamVL = yopt(2*n_units+n_variable+1:2*(n_units+n_variable));
if n_extra ~= 0
   lamExtra = yopt(2*(n_units+n_variable)+1:end-1);
end
maxSensitivity.expu = zeros(n_units,1);
maxSensitivity.expl = zeros(n_units,1);
for i1 = 1:n_units
   if isa(obj.DatasetUnits.Values(i1).SurrogateModel,'B2BDC.B2Bmodels.QModel')
      maxSensitivity.expu(i1) = lamEU(i1)*d(i1)*S(1,1);
      maxSensitivity.expl(i1) = lamEL(i1)*d(i1)*S(1,1);
   elseif isa(obj.DatasetUnits.Values(i1).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
      tmpD = obj.DatasetUnits.Values(i1).SurrogateModel.NormalizedDenominator;
      maxSensitivity.expu(i1) = lamEU(i1)*d(i1)*trace(tmpD*S);
      maxSensitivity.expl(i1) = lamEL(i1)*d(i1)*trace(tmpD*S);
   end
end
maxSensitivity.varu = -1.05*lamVU.*vd.^2;
maxSensitivity.varl = -1.05*lamVL.*vd.^2;
if n_extra ~= 0
   count = 1;
   for k1 = 1:n_extra
      i1 = extraIdx(k1,1);
      j1 = extraIdx(k1,2);
      tmpLam = lamExtra(count:count+3);
      maxSensitivity.varu(i1) = maxSensitivity.varu(i1)-vd(i1)*(tmpLam(3)*(vl(j1)*S(1,1)-S(1,j1))+tmpLam(4)*(-vu(j1)*S(1,1)+S(1,j1)));
      maxSensitivity.varu(j1) = maxSensitivity.varu(j1)-vd(j1)*(tmpLam(2)*(vl(i1)*S(1,1)-S(1,i1))+tmpLam(4)*(-vu(i1)*S(1,1)+S(1,i1)));
      maxSensitivity.varl(i1) = maxSensitivity.varl(i1)-vd(i1)*(tmpLam(1)*(vl(j1)*S(1,1)-S(1,j1))+tmpLam(2)*(-vu(j1)*S(1,1)+S(1,j1)));
      maxSensitivity.varl(j1) = maxSensitivity.varl(j1)-vd(j1)*(tmpLam(1)*(vl(i1)*S(1,1)-S(1,i1))+tmpLam(3)*(-vu(i1)*S(1,1)+S(1,i1)));
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
         targetMatrix = -targetN;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         Coef0 = model0.NormalizedCoefMatrix;
         targetMatrix = zeros(n_variable+1,n_variable+1);
         for j = 1:length(id2)
            for k = 1:length(id2)
               targetMatrix(id1(j),id1(k)) = -Coef0(id2(j),id2(k));
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
         As(:,end) = -targetD;
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         As(1,end) = -1;
      end
      At = [Alam; As];
      c = [clam; cs];
      b = spalloc(n_opt,1,1);
      b(end) = -1;
      K.l = n_opt-1;
      K.s = n_variable+1;
      A = At';    
   end

end