function y = sedumiquadfitinfnorm(X,Y,vars,nCV)
% Quadratic fit of the given data X(nSample-by-nVariable) and
% Y(nSample-by-1) with variable information in the
% B2BDC.B2Bvariables.VariableList object vars. The fitting criteria is to
% minimize the infinity norm of error
% output:
% A B2BDC.B2Bmodels.QModel object

% Created: July 22, 2015    Wenyu Li


[n_sample, n_variable] = size(X);
if ~isvector(Y) || n_sample ~= length(Y)
   error('Wrong dimension of input data')
end
my = mean(Y);
dy = 0.5*(max(Y)-min(Y));
if nCV == 1
   % mx = mean(vars.calBound')';
   % dx = 0.5*diff(vars.calBound')';
   % Tx = [1, zeros(1,n_variable);
   %    -mx./dx , diag(1./dx)];
   % Y = (Y-my)/dy;
   % x1 = [ones(n_sample,1),X];
   % x2 = Tx*x1';
   % x2 = x2';
   % Xn = x2(:,2:end);
   % xNew = B2BDC.Fitting.expandBasis(Xn);
   xNew = B2BDC.Fitting.expandBasis(X);
   n_Coef = size(xNew,2);
   n_opt = n_Coef+1+2*n_sample;
   A = zeros(2*n_sample, n_opt);
   A(1:n_sample,n_Coef+1) = -ones(n_sample,1);
   A(n_sample+1:end,n_Coef+1) = ones(n_sample,1);
   A(1:n_sample,n_Coef+2:n_Coef+1+n_sample) = eye(n_sample);
   A(n_sample+1:end,n_Coef+2+n_sample:end) = -eye(n_sample);
   A(:,1:n_Coef) = repmat(xNew,2,1);
   b = repmat(Y,2,1);
   c = zeros(n_opt,1);
   c(n_Coef+1) = 1;
   k.f = n_Coef;
   k.l = 2*n_sample+1;
   pars.fid = 0;
   [x,y,info] = sedumi(A,b,c,k,pars);
   coefVec = x(1:n_Coef);
   c1 = B2BDC.Fitting.vec2coef(coefVec,n_variable);
   yScale.my = my;
   yScale.dy = dy;
   y = B2BDC.B2Bmodels.QModel(c1,vars,yScale);
   yhat = xNew * coefVec;
   err = abs(yhat-Y);
   y.ErrorStats.absMax = max(err);
   y.ErrorStats.absAvg = mean(err);
   err = err./abs(Y);
   y.ErrorStats.relMax = max(err);
   y.ErrorStats.relAvg = mean(err);
else
   nTest = floor(n_sample/nCV);
   err.absMax = 0;
   err.absAvg = 0;
   err.relMax = 0;
   err.relAvg = 0;
   for i = 1:nCV
      if i ~= nCV
         idTest = (1:nTest)+(i-1)*nTest;
      else
         idTest = (nCV-1)*nTest+1:n_sample;
      end
      idTrain = setdiff(1:n_sample,idTest);
      xTrain = X(idTrain,:);
      yTrain = Y(idTrain);
      n1 = length(yTrain);
      xTest = X(idTest,:);
      yTest = Y(idTest);
      xFit = B2BDC.Fitting.expandBasis(xTrain);
      xTest = B2BDC.Fitting.expandBasis(xTest);
      n_Coef = size(xFit,2);
      n_opt = n_Coef+1+2*n1;
      A = zeros(2*n1, n_opt);
      A(1:n1,n_Coef+1) = -ones(n1,1);
      A(n1+1:end,n_Coef+1) = ones(n1,1);
      A(1:n1,n_Coef+2:n_Coef+1+n1) = eye(n1);
      A(n1+1:end,n_Coef+2+n1:end) = -eye(n1);
      A(:,1:n_Coef) = repmat(xFit,2,1);
      b = repmat(yTrain,2,1);
      c = zeros(n_opt,1);
      c(n_Coef+1) = 1;
      k.f = n_Coef;
      k.l = 2*n1+1;
      pars.fid = 0;
      [x,y,info] = sedumi(A,b,c,k,pars);
      coefVec = x(1:n_Coef);
      dy = abs(xTest*coefVec-yTest);
      err.absMax = max(err.absMax,max(dy));
      err.absAvg = err.absAvg+mean(dy);
      dy = dy./abs(yTest);
      err.relMax = max(err.relMax,max(dy));
      err.relAvg = err.relAvg+mean(dy);
   end
   xNew = B2BDC.Fitting.expandBasis(X);
   n_Coef = size(xNew,2);
   n_opt = n_Coef+1+2*n_sample;
   A = zeros(2*n_sample, n_opt);
   A(1:n_sample,n_Coef+1) = -ones(n_sample,1);
   A(n_sample+1:end,n_Coef+1) = ones(n_sample,1);
   A(1:n_sample,n_Coef+2:n_Coef+1+n_sample) = eye(n_sample);
   A(n_sample+1:end,n_Coef+2+n_sample:end) = -eye(n_sample);
   A(:,1:n_Coef) = repmat(xNew,2,1);
   b = repmat(Y,2,1);
   c = zeros(n_opt,1);
   c(n_Coef+1) = 1;
   k.f = n_Coef;
   k.l = 2*n_sample+1;
   pars.fid = 0;
   [x,y,info] = sedumi(A,b,c,k,pars);
   coefVec = x(1:n_Coef);
   c1 = B2BDC.Fitting.vec2coef(coefVec,n_variable);
   yScale.my = my;
   yScale.dy = dy;
   y = B2BDC.B2Bmodels.QModel(c1,vars,yScale);
   y.ErrorStats.absMax = err.absMax;
   y.ErrorStats.absAvg = err.absAvg/nCV;
   y.ErrorStats.relMax = err.relMax;
   y.ErrorStats.relAvg = err.relAvg/nCV;
end