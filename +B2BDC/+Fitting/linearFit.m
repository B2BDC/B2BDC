<<<<<<< Updated upstream
function y = linearFit(X,Y,vars)
% Linear fit of the given data X(nSample-by-nVariable) and
% Y(nSample-by-1) with variable information in the
% B2BDC.B2Bvariables.VariableList object vars. The fitting criteria is to
% minimize the 2-norm of error
% output:
% A B2BDC.B2Bmodels.QModel object

xNew = [ones(size(X,1),1) X];
yScale.my = mean(Y);
yScale.dy = 0.5*(max(Y)-min(Y));
coefVec = xNew \ Y;
y = B2BDC.B2Bmodels.LModel(coefVec, vars, yScale);
yhat = xNew*coefVec;
err = abs(yhat-Y);
y.ErrorStats.absMax = max(err);
y.ErrorStats.absAvg = mean(err);
err = err./abs(Y);
y.ErrorStats.relMax = max(err);
=======
function y = linearFit(X,Y,vars)
% Linear fit of the given data X(nSample-by-nVariable) and
% Y(nSample-by-1) with variable information in the
% B2BDC.B2Bvariables.VariableList object vars. The fitting criteria is to
% minimize the 2-norm of error
% output:
% A B2BDC.B2Bmodels.QModel object

xNew = [ones(size(X,1),1) X];
yScale.my = mean(Y);
yScale.dy = 0.5*(max(Y)-min(Y));
coefVec = xNew \ Y;
y = B2BDC.B2Bmodels.LModel(coefVec, vars, yScale);
yhat = xNew*coefVec;
err = abs(yhat-Y);
y.ErrorStats.absMax = max(err);
y.ErrorStats.absAvg = mean(err);
err = err./abs(Y);
y.ErrorStats.relMax = max(err);
>>>>>>> Stashed changes
y.ErrorStats.relAvg = mean(err);