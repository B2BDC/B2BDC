function y = sedumiquadfit2norm(X,Y,vars)
% Quadratic fit of the given data X(nSample-by-nVariable) and
% Y(nSample-by-1) with variable information in the
% B2BDC.B2Bvariables.VariableList object vars. The fitting criteria is to
% minimize the 2-norm of error
% output:
% A B2BDC.B2Bmodels.QModel object

% Created: July 22, 2015    Wenyu Li

[n_sample, n_variable] = size(X);
mx = mean(vars.calBound')';
dx = 0.5*diff(vars.calBound')';
Tx = [1, zeros(1,n_variable);
   -mx./dx , diag(1./dx)];
if ~isvector(Y) || n_sample ~= length(Y)
   error('Wrong dimension of input data')
end
my = mean(Y);
dy = 0.5*(max(Y)-min(Y));
Y = (Y-my)/dy;
x1 = [ones(n_sample,1),X];
x2 = Tx*x1';
x2 = x2';
Xn = x2(:,2:end);
xNew = B2BDC.Fitting.expandBasis(Xn);
n_Coef = size(xNew,2);
n_opt = n_Coef+1+n_sample;
A = zeros(n_sample,n_opt);
b = Y;
A(:,1:n_Coef) = xNew;
A(:,n_Coef+2:end) = -eye(n_sample);
c = zeros(n_opt,1);
c(n_Coef+1) = 1;
k.f = n_Coef;
k.q = n_sample+1;
pars.fid = 0;
[x,y,info] = sedumi(A,b,c,k,pars);
coefVec = x(1:n_Coef);
c1 = B2BDC.Fitting.vec2coef(coefVec,n_variable);
yScale.my = my;
yScale.dy = dy;
y = B2BDC.B2Bmodels.QModel(c1,vars,yScale);
yhat = xNew * coefVec;
err = dy*abs(yhat-Y);
y.ErrorStats.absMax = max(err);
y.ErrorStats.absAvg = mean(err);
err = err./abs(Y*dy+my);
y.ErrorStats.relMax = max(err);
y.ErrorStats.relAvg = mean(err);
