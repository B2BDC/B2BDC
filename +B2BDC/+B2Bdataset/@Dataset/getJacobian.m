function J = getJacobian(obj)
% Calculate the indication n_variable-by-n_variable matrix J of a
% B2BDC.B2Bdataset.Dataset object. n_variable is the total number of
% variables in the dataset. J is in upper-triangular form so J_ij is an
% indication about the total influence of the interaction between x_i and
% x_j on the dataset.

% Created: July 14, 2015     Wenyu Li
%  Modified: August 5, 2015    Wenyu Li  (Works for both RQModel and QModel)

n_units = obj.Length;
units = obj.DatasetUnits.Values;
allVarname = obj.VarNames;
n_variable = length(allVarname);
J = zeros(n_variable,n_variable);
for i = 1:n_units
   unitModel = units(i).SurrogateModel;
   [~,id1,id2] = intersect(allVarname,unitModel.VarNames);
   if isa(unitModel,'B2BDC.B2Bmodels.RQModel')
      NumCoef0 = zeros(n_variable,n_variable);
      NumCoef1 = 2*unitModel.NormalizedNumerator(2:end,2:end);
      DenCoef0 = zeros(n_variable,n_variable);
      DenCoef1 = 2*unitModel.NormalizedDenominator(2:end,2:end);
      for j = 1:length(id2)
         for k = 1:length(id2)
            NumCoef0(id1(j),id1(k)) = NumCoef1(id2(j),id2(k));
            DenCoef0(id1(j),id1(k)) = DenCoef1(id2(j),id2(k));
         end
      end
      J = J + abs(NumCoef0) + abs(DenCoef0);
   elseif isa(unitModel,'B2BDC.B2Bmodels.QModel')
      quadCoef0 = zeros(n_variable,n_variable);
      quadCoef1 = unitModel.Hessian;
      for j = 1:length(id2)
         for k = 1:length(id2)
            quadCoef0(id1(j),id1(k)) = quadCoef1(id2(j),id2(k));
         end
      end
      J = J + abs(quadCoef0);
   end
end
for i = 1:n_variable
   J(i:end,i) = zeros(n_variable-i+1,1);
end