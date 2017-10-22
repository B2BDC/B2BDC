<<<<<<< Updated upstream
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
=======
function xopt = optimizationParameter(obj,opt,id,wY,lambda,wX,alpha)
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
if nargin < 4 || isempty(wY)
   switch opt.ConsistencyMeasure
      case 'relative'
         wY = 1./obj.calObserve;
      case 'absolute'
         wY = ones(n_units,1);
      otherwise
         qbd = obj.calBound;
         wY = 1./(qbd(:,2)-qbd(:,1));
   end
else
   if length(wY) ~= n_units
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
if nargin < 6 || isempty(wX)
   switch opt.ConsistencyMeasure
      case 'relative'
         xbd = obj.Variables.calBound;
         wX = 1./(xbd(:,2)-xbd(:,1));
      otherwise
         wX = ones(nVar,1); 
   end 
else
   if length(wX) ~= nVar
      error('Invalid length of variable weights')
   end
end
if nargin < 7
   alpha = 0;
end
% flag = quadratictest(obj);
if Nflag
   xopt = optParameter_fmincon(obj,id,Nflag,Hflag,lambda,opt);
   % xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,lambda,opt);
elseif any(wX)
   xopt = optParameter_fmincon_withPenalty(obj,id,Hflag,wX,wY,alpha,opt);
else
   xopt = optParameter_fmincon(obj,id,Nflag,Hflag,wY,opt);
   % xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,w,opt);
end
>>>>>>> Stashed changes
