function flag = quadratictest(obj)
% Test if all the surrogate models in the dataset are quadratic models.
% Return true if they are, false otherwise.

% Created June 30, 2015     Wenyu Li

flag = true;
dsUnits = obj.DatasetUnits.Values;
for i = 1:length(dsUnits)
   testModel = dsUnits(i).SurrogateModel;
   if ~isa(testModel,'B2BDC.B2Bmodels.QModel')
      flag = false;
      break;
   end
end
end