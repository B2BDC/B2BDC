function rqModel = sedumirqfitinfnorm(X,Y,vars,b2bopt,K)
% fit the data with a rational quadratic function. X is a
% nSample-by-nVariable matrix and Y is a nSample-by-1 column vector. vars
% is a B2BDC.B2Bvariables.VariableList object contains information of
% model variables. K is an optional input, which is the upper
% bound of the denominator over its domain and must be larger than one. If
% K is not provided, the function itself will implement a 4-fold cross
% validation to select a best K value between 2 and 20. It returns a
% B2BDC.B2Bmodel.RQModel object
% Input format for sedumi solver is created based on the dual variable 
%  y = [ n  D  lambda  gamma r]' :
%  N - Vector coefficient for numerator quadratics
%  D - Vector coefficient for denominator quadratics
% lambda - Laguagian multipliers for box constraint of x with D > 1
% gamma  - Laguagian multipliers for box constraint of x with D < K
%  r - A variable used to create the semidefinite constraint


% output:
% A B2BDC.B2Bmodels.RQModel object

% Created: July 22, 2015    Wenyu Li
%  Modified: August 16, 2015  Wenyu Li

[n_sample,n_variable] = size(X);
% mx = mean(vars.calBound')';
% dx = 0.5*diff(vars.calBound')';
% Tx = [1, zeros(1,n_variable);
%    -mx./dx , diag(1./dx)];
if ~isvector(Y) || n_sample ~= length(Y)
   error('Wrong dimension of input data')
end
my = mean(Y);
dy = 0.5*(max(Y)-min(Y));
yScale.my = my;
yScale.dy = dy;
% Y = (Y-my)/dy;
% x1 = [ones(n_sample,1) X];
% x2 = Tx*x1';
% x2 = x2';
X0 = X;
% X = x2(:,2:end);
nc = 0.5*(n_variable+1)*(n_variable+2);
LB = -ones(n_variable,1);
UB = -LB;
if nargin == 4
   if b2bopt.Display
      disp('Start cross validation to select best K value...')
   end
   if 0.75*n_sample > 2*nc
      n_cv = 4;
   elseif 0.9*n_sample > 2*nc
      n_cv = 10;
   else
      error('Too less sample point to fit a rational quadratic model')
   end   
   n_validate = floor(n_sample/n_cv);
   n_k = 7;
   K_option = linspace(2,20,n_k);
   Nvec = cell(n_k,n_cv);
   Dvec = cell(n_k,n_cv);
   result = zeros(n_k,n_cv);
   for j1 = 1:n_k
      k1 = K_option(j1);
      if b2bopt.Display
         disp(['Trying K value ' int2str(k1) ', ' int2str(n_k-j1) ' left'])
      end
      for j2 = 1:n_cv
         idx = true(n_sample,1);
         idx((j2-1)*n_validate+1:j2*n_validate) = false;
         xtrain = X(idx,:);
         xtest = X(~idx,:);
         ytrain = Y(idx);
         ytest = Y(~idx);
         t_1 = 0;
         yy = B2BDC.Fitting.sedumiquadfitinfnorm(X0(idx,:),ytrain,vars);
         t_2 = 1.1*yy.ErrorStats.absMax;
         tolerance = 1e-4*t_2;
         while t_2-t_1 > tolerance
            t_m = 0.5*(t_1+t_2);
            [A,b,c,k] = setSedumiCoef(xtrain,ytrain,k1);
            pars.fid = 0;
            pars.eps = 1e-7;
            [~,yopt] = sedumi(A,b,c,k,pars);
            if yopt(end) >= 0
               t_2 = t_m;
               Nvec{j1,j2} = yopt(1:nc);
               Dvec{j1,j2} = yopt(nc+1:2*nc);
            else
               t_1 = t_m;
            end
         end
         x = B2BDC.Fitting.expandBasis(xtest);
         if ~isempty(Nvec{j1,j2})
            fx = x*Nvec{j1,j2}./(x*Dvec{j1,j2})-ytest;
         else
            fx = ytest;
         end
         result(j1,j2) = max(abs(fx));
      end
   end
   if b2bopt.Display
      disp('Cross validation finished')
   end
   result_cv = mean(result,2);
   [~,idx] = min(result_cv);
   kbest = K_option(idx);
   kk = kbest;
elseif nargin == 5 && K > 1
   if n_sample > 2*nc
      kk = K;
   else
      error('Too less sample point to fit a rational quadratic model')
   end
end
if b2bopt.Display
   disp(['Best K value is ' int2str(kk) ', start fitting with the best K value'])
end
t_1 = 0;
yy = B2BDC.Fitting.sedumiquadfitinfnorm(X0,Y,vars);
t_2 = 1.1*yy.ErrorStats.absMax;
tolerance = 1e-4*t_2;
while t_2-t_1>tolerance
   t_m = 0.5*(t_1+t_2);
   [A,b,c,k] = setSedumiCoef(X,Y,kk);
   pars.fid = 0;
   pars.eps = 1e-7;
   [~,yopt] = sedumi(A,b,c,k,pars);
   if yopt(end) >= 0
      t_2 = t_m;
      Nvec = yopt(1:nc);
      Dvec = yopt(nc+1:2*nc);
   else
      t_1 = t_m;
   end
end
Nmod = B2BDC.Fitting.vec2coef(Nvec,n_variable);
Dmod = B2BDC.Fitting.vec2coef(Dvec,n_variable);
rqModel = B2BDC.B2Bmodels.RQModel(Nmod,Dmod,vars,kk,yScale);
x = B2BDC.Fitting.expandBasis(X);
err = abs(x*Nvec./(x*Dvec)-Y);
rqModel.ErrorStats.absMax = max(err);
rqModel.ErrorStats.absAvg = mean(err);
err = err./abs(Y);
rqModel.ErrorStats.relMax = max(err);
rqModel.ErrorStats.relAvg = mean(err);
if b2bopt.Display
   disp('Done');
end


   
   
   
   
   function [A,b,c,k] = setSedumiCoef(Xsdm,Ysdm,Kval)
      xdata = B2BDC.Fitting.expandBasis(Xsdm);
      n_opt = 2*(nc+n_variable)+1;
      n_data = size(xdata,1);
      Ac1 = [xdata, repmat(-(Ysdm+t_m),1,nc).*xdata, zeros(n_data,2*n_variable+1)];
      cc1 = spalloc(n_data,1,0);
      Ac2 = [-xdata, repmat(Ysdm-t_m,1,nc).*xdata, zeros(n_data,2*n_variable+1)];
      cc2 = spalloc(n_data,1,0);
      Dmat = zeros((n_variable+1)^2,nc);
      Dmat(1,1) = 1;
      for i = 1:n_variable
         tep = zeros(n_variable+1);
         tep(1,i+1) = 0.5;
         tep(i+1,1) = 0.5;
         Dmat(:,i+1) = vec(tep);
      end
      count = n_variable+2;
      for i1 = 2:n_variable+1
         for i2 = i1:n_variable+1
            tep = zeros(n_variable+1);
            tep(i1,i2) = 0.5;
            Dmat(:,count) = vec(tep+tep');
            count = count+1;
         end
      end
      Hmat = zeros((n_variable+1)^2,n_variable);
      for i = 1:n_variable
         tep = zeros(n_variable+1);
         tep(1,1) = LB(i)*UB(i);
         tep(1,i+1) = -0.5*(LB(i)+UB(i));
         tep(i+1,1) = -0.5*(LB(i)+UB(i));
         tep(i+1,i+1) = 1;
         Hmat(:,i) = vec(tep);
      end
      As1 = [spalloc((n_variable+1)^2,nc,0), -Dmat, -Hmat, spalloc((n_variable+1)^2,n_variable,0), vec(speye(n_variable+1))];
      cs1 = zeros((n_variable+1)^2,1);
      cs1(1) = -1;
      As2 = [spalloc((n_variable+1)^2,nc,0), Dmat, spalloc((n_variable+1)^2,n_variable,0), -Hmat, vec(speye(n_variable+1))];
      cs2 = zeros((n_variable+1)^2,1);
      cs2(1) = Kval;
      Alam = [spalloc(2*n_variable,2*nc,0), -speye(2*n_variable), spalloc(2*n_variable,1,0)];
      clam = spalloc(2*n_variable,1,0);
      b = spalloc(n_opt,1,0);
      b(end) = 1;
      At = [Ac1; Ac2; Alam; As1; As2];
      c = [cc1; cc2; clam; cs1; cs2];
      A = At';
      k.l = 2*(n_data+n_variable);
      k.s = [n_variable+1, n_variable+1];
   end
end