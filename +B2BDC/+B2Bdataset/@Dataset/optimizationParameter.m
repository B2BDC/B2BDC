function xopt = optimizationParameter(obj,opt,id,w,lambda)
%XOPT = OPTIMIZATIONONF(OBJ,OPT,W,LAMBDA) returns a set of optimal variables with
%certain optimization criteria specified in opt.

% Created: July 19, 2016

if nargin < 2
   opt = generateOpt('Display',false);
end
if nargin < 3
   id = [];
end
Nflag = strcmp(opt.Optimization,'1N-F');
Hflag = strcmp(opt.Optimization,'LS-H');
if ~Hflag
   if ~obj.isConsistent(opt)
      error('The dataset is not consistent');
   end
end
n_units = obj.Length;
nVar = obj.Variables.Length;
if nargin < 4 || isempty(w)
   switch opt.ConsistencyMeasure
      case 'relative'
         w = 1./obj.calObserve;
      case 'absolute'
         w = ones(n_units,1);
      otherwise
         qbd = obj.calBound;
         w = 1./(qbd(:,2)-qbd(:,1));
   end
else
   if length(w) ~= n_units
      error('Invalid length of QOI weights')
   end
end
if nargin < 5 || isempty(lambda)
   switch opt.ConsistencyMeasure
      case 'relative'
         xbd = obj.Variables.calBound;
         lambda = 1./(xbd(:,2)-xbd(:,1));
      otherwise
         lambda = ones(nVar,1);    
   end 
else
   if length(lambda) ~= nVar
      error('Invalid length of variable weights')
   end
end
flag = quadratictest(obj);
if Nflag
   xopt = optParameter_fmincon(obj,id,flag,Nflag,Hflag,lambda,opt);
   % xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,lambda,opt);
else
   xopt = optParameter_fmincon(obj,id,flag,Nflag,Hflag,w,opt);
   % xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,w,opt);
end
