function newPoly = expandDimension(obj,newVar)
% Returns a new polynomial model with its basis's dimension expanded to a
% new variable set specified by the input new VariableList.

%  Created: Nov 30, 2015     Wenyu Li

oldBasis = obj.Basis;
newBasis = oldBasis.expandDimension(obj.Variables,newVar);
Coef = obj.Coefficient;
dataStat = obj.DataStats;
err = obj.ErrorStats;
newPoly = B2BDC.B2Bmodels.PolyModel(newBasis, Coef, newVar, dataStat, err);