function y = generateModelbyFit(X,Y,vars,mtype,opt,kk)
% Returns a B2BDC.B2Bmodels.Model object by fitting the data X, Y
% with possibly different choice of model types. The last 4 input are
% optional. If they are not empty, the default for norflagx and norflagy  
% are true and the default for mtype is 'qinf'.
% The input arguments are:
% 
%   X    -  nSample-by-nVar matrix of x data
%   Y    -  nSample-by-1 colum vector of y data
%  vars  -  A B2BDC.B2Bvariables.VariableList object contains variable
%           informations (name, domain etc.)
%  mtype -  The model type that fits the data defined as following strings
%           'qinf': qudaratic model with infinity norm error minimized
%           'q2norm': quadratic model with 2 norm error minimized
%           'rq': rational quadratic model with infinity norm error minimized
%   kk   -  The parameter used to upper bound the denominator of ratioanl
%           quadratic model over its variable domain 

%  Created: Nov 1, 2015    Wenyu Li

if nargin < 3
   error('Not enough input variables')
end
if nargin > 3 && ~isempty(mtype)
   switch mtype
      case 'qinf'
         y = B2BDC.Fitting.sedumiquadfitinfnorm(X,Y,vars);
      case 'q2norm'
         y = B2BDC.Fitting.sedumiquadfit2norm(X,Y,vars);
      case 'rq'
         if nargin < 5
            opt = generateOpt('Display',false);
         end
         if nargin > 5 && ~isempty(kk)
            y = B2BDC.Fitting.sedumirqfitinfnorm(X,Y,vars,opt,kk);
         else
            y = B2BDC.Fitting.sedumirqfitinfnorm(X,Y,vars,opt);
         end
   end
else
   y = B2BDC.Fitting.sedumiquadfitinfnorm(X,Y,vars);
end