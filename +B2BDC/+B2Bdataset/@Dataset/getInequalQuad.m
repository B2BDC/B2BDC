function [Qunits, Qx, Qextra, n_extra, extraIdx] = getInequalQuad(obj,bds,frac)
% Return the quadratic form of inequality constraints of the
% B2BDC.B2Bdataset.Dataset object. The returned quadratic form matrix is
% with respect to all the active variables of the dataset. 
% Input:
%  bds  -  A nUnits-by-2 matrix defines the lower and upper bounds for each
%          surrogate model in the dataset unit.
%  frac -  Fraction of extra constraints used in the optimization
%          0 < frac < 100, if frac == -1, then automatically linear paris
%          with influence factor greater than 5% of the most influential pair
%          will be included
% Output:
%  Qunits - Inequality quadratic form of the observation QOI to be within
%           the uncertainty bounds
%    Qx   - Inequality quadratic form of the linear box constraints of all
%           the active variables of the dataset
%  Qextra - Inequality quadratic form of extra combination pairs of linear
%           constraints on the active variables of the dataset
%  n_extra - Number of extra variable pairs used in the SDP
% extraIdx - A n_extra-by-2 matrix stores the index of each variable pair
% All the inequality quadratic matrix is formed as Q <= 0

% Created: August 10, 2015     Wenyu Li

vars = obj.Variables.Values;
varName = obj.VarNames;
n_variable = length(vars);
LB = -ones(n_variable,1);
UB = -LB;
units = obj.DatasetUnits.Values;
n_units = length(units);
J = obj.getJacobian;
[Jsort,idx] = sort(J(:),'descend');
if frac == -1
   n_extra = sum(Jsort(:) >= 0.05*Jsort(1));
else
   n_extra = floor(0.005*n_variable*(n_variable-1)*frac);
end
Qunits = cell(2*n_units,1);
Qx = cell(2*n_variable,1);
Qextra = cell(4*n_extra,1);
for i = 1:n_units
   model = units(i).SurrogateModel;
   modelVar = model.VarNames;
   [~,id1,id2] = intersect(varName,modelVar);
   id1 = [1;id1+1];
   id2 = [1;id2+1];
   if isa(model,'B2BDC.B2Bmodels.RQModel')
      N = model.Numerator;
      D = model.Denominator;
      constrN = zeros(n_variable+1,n_variable+1);
      constrD = zeros(n_variable+1,n_variable+1);
      for j = 1:length(id2)
         for k = 1:length(id2)
            constrN(id1(j),id1(k)) = N(id2(j),id2(k));
            constrD(id1(j),id1(k)) = D(id2(j),id2(k));
         end
      end
      Qunits{2*i-1} = constrN-bds(i,2)*constrD;
      Qunits{2*i} = bds(i,1)*constrD-constrN;
   elseif isa(model,'B2BDC.B2Bmodels.QModel')
      Coef = model.CoefMatrix;
      constrMatrix1 = zeros(n_variable+1,n_variable+1);
      constrMatrix2 = zeros(n_variable+1,n_variable+1);
      for j = 1:length(id2)
         for k = 1:length(id2)
            constrMatrix1(id1(j),id1(k)) = Coef(id2(j),id2(k));
            constrMatrix2(id1(j),id1(k)) = -Coef(id2(j),id2(k));
         end
      end
      constrMatrix1(1,1) = constrMatrix1(1,1)-bds(i,2);
      constrMatrix2(1,1) = constrMatrix2(1,1)+bds(i,1);
      Qunits{2*i-1} = constrMatrix1;
      Qunits{2*i} = constrMatrix2;
   end
end
eta = 0.1;
for i = 1:n_variable
   lx = LB(i);
   ux = UB(i);
   Qx{2*i-1} = sparse([1,1,i+1,i+1],[1,i+1,1,i+1],...
      [(lx-eta)*ux,-(lx+ux-eta)/2,-(lx+ux-eta)/2,1],n_variable+1,n_variable+1);
   Qx{2*i} = sparse([1,1,i+1,i+1],[1,i+1,1,i+1],...
      [lx*(ux+eta),-(lx+ux+eta)/2,-(lx+ux+eta)/2,1],n_variable+1,n_variable+1);
end
extraIdx = zeros(n_extra,2);
if n_extra ~= 0
   count = 1;
   for k = 1:n_extra
      [i,j] = ind2sub(size(J),idx(k));
      extraIdx(k,:) = [i,j];
      lxi = LB(i);
      uxi = UB(i);
      lxj = LB(j);
      uxj = UB(j);
      Qextra{count} = sparse([1,1,i+1,1,j+1,i+1,j+1],[1,i+1,1,j+1,1,j+1,i+1],...
         [-lxi*lxj,0.5*lxj,0.5*lxj,0.5*lxi,0.5*lxi,-0.5,-0.5],n_variable+1,n_variable+1);
      Qextra{count+1} = sparse([1,1,i+1,1,j+1,i+1,j+1],[1,i+1,1,j+1,1,j+1,i+1],...
         [lxi*uxj,-0.5*uxj,-0.5*uxj,-0.5*lxi,-0.5*lxi,0.5,0.5],n_variable+1,n_variable+1);
      Qextra{count+2} = sparse([1,1,i+1,1,j+1,i+1,j+1],[1,i+1,1,j+1,1,j+1,i+1],...
         [uxi*lxj,-0.5*lxj,-0.5*lxj,-0.5*uxi,-0.5*uxi,0.5,0.5],n_variable+1,n_variable+1);
      Qextra{count+3} = sparse([1,1,i+1,1,j+1,i+1,j+1],[1,i+1,1,j+1,1,j+1,i+1],...
         [-uxi*uxj,0.5*uxj,0.5*uxj,0.5*uxi,0.5*uxi,-0.5,-0.5],n_variable+1,n_variable+1);
      count = count+4;
   end
end