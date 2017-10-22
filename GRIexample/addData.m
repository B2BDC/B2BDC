function  y = addData(b2bDS,expData,modelData)
% Create dataset unit from experimental data and model data and then add to
% the B2BDC dataset

% Created: June 23, 2015    Wenyu Li

n_units = size(expData,1);
for i = 1:n_units
   modelVar = modelData{i,3};
   n = size(modelVar,1);
   varList = B2BDC.B2Bvariables.VariableList();
   varLB = zeros(n,1);
   varUB = zeros(n,1);
   for j = 1:n
      modelVariable = B2BDC.B2Bvariables.ModelVariable(modelVar{j,1},...
         modelVar{j,2},modelVar{j,3},modelVar{j,4});
      varList = varList.add(modelVariable);
      varLB(j) = modelVar{j,2};
      varUB(j) = modelVar{j,3};
   end
   modelM = modelData{i,2};
   switch modelM{1}
      case 'quadratic'
         if length(modelM) == 2
            quadCoef = modelM{2};
            unitModel = B2BDC.B2Bmodels.QModel(quadCoef,varList);
         elseif length(modelM) == 3
            Xdata = modelM{2};
            Ydata = modelM{3};
            unitModel = B2BDC.Fitting.sedumiquadfitinfnorm(Xdata,Ydata,[varLB,varUB]);
         end
      case 'rational quadratic'
         if length(modelM) == 4
            matrixNum = modelM{2};
            matrixDen = modelM{3};
            kVal = modelM{4};
            unitModel = B2BDC.B2Bmodels.RQModel(matrixNum,matrixDen,kVal,varList);
         elseif length(modelM) == 3
            Xdata = modelM{2};
            Ydata = modelM{3};
            unitModel = B2BDC.Fitting.sedumirqfitinfnorm(Xdata,Ydata,varList,kVal);
         end
   end
   unitName = expData{i,1};
   unitLB = expData{i,2};
   unitUB = expData{i,3};
   if size(expData,2) == 3
      unitVal = 0.5*(expData{i,2}+expData{i,3});
   elseif size(expData,2) == 4
      unitVal = expData{i,4};
   end
   dsUnit = B2BDC.B2Bdataset.DatasetUnit(unitName,unitModel,unitLB,unitUB,unitVal);
   b2bDS.addDSunit(dsUnit);
end
y = b2bDS;