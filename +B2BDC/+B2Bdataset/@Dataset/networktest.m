function flag = networktest(obj)
% Test if all the surrogate models in the dataset are network models.
% Return true if they are, false otherwise.

% Created Aug 21, 2018     Wenyu Li

flag = true;
dsUnits = obj.DatasetUnits.Values;
for i = 1:length(dsUnits)
   testModel = dsUnits(i).SurrogateModel;
   if ~isa(testModel,'B2BDC.B2Bmodels.NNModel')
      flag = false;
      break;
   end
end