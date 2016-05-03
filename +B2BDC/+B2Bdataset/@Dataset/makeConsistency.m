function makeConsistency(obj,itMax,opt)
% Updating an inconsistent dataset iteratively into consistent by gradually
% expanding the experiment bounds selected by ranked normalized
% sensitivity. In each iteration step, the selected dataset unit bounds are
% expanded by 5% of the current uncertainty (difference between upper and
% lower bounds). 
% The input arguments are:
%    opt - B2BDC.Option object
%    itMax - Maximum iteration number (optional). If not specified, the
%            default value is 10.

% Created: August 25, 2015    Wenyu Li

if nargin < 2
   itMax = 10;
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
end
iter = 1;
deta = 0;
n_Unit = obj.Length;
while ~obj.isConsistent(opt)
   c = obj.ConsistencyMeasure;
   display(['Consistency measure is: [' num2str(c(1)) ', ' num2str(c(2)) ']'])
   display(['Iteration count: ' int2str(iter)]);
   bds = obj.calBound;
   lb = bds(:,1);
   ub = bds(:,2);
   db = diff(bds');
   ob = obj.calObserve;
   s = [obj.ConsistencySensitivity.expu; obj.ConsistencySensitivity.expl];
   [~,id] = sort(s,'descend');
   if id(1) <= n_Unit
      dsName = [obj.DatasetUnits.Values(id(1)).Name ' UB'];
   else
      dsName = [obj.DatasetUnits.Values(id(1)-n_Unit).Name ' LB'];
   end
   display(['Change the dataset unit ' dsName]);
   deta = deta + 0.05*db(id(1));
   newd = db(id(1))*1.05;
   if id(1) <= n_Unit
      newbd = [lb(id(1)), lb(id(1))+newd, ob(id(1))];
   else
      newbd = [ub(id(1)-n_Unit)-newd, ub(id(1)-n_Unit), ob(id(1)-n_Unit)];
   end
   if id(1) <= n_Unit
      obj.changeBounds(newbd,id(1));
   else
      obj.changeBounds(newbd,id(1)-n_Unit);
   end
   iter = iter+1;
   if iter >= itMax
      break
   end
end
disp('done');
if obj.isConsistent(opt)
   c = obj.ConsistencyMeasure;
   display(['Consistency measure is: [' num2str(c(1)) ', ' num2str(c(2)) ']'])
   disp(['The total adjustment in experiment bounds are: ' num2str(deta)])
else
   disp(['The dataset fail to get consistent after ' num2str(itMax) ' iterations'])
end