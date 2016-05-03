function evalConsistencyrel(obj,b2bopt)
% Calculate the inner and outer bound for the consistency measure (in 
% percentage) of a B2BDC.B2Bdataset.Dataset object
% frac - Fraction of extra constraints used in the optimization
%        0 < frac < 100, if frac == -1, then automatically linear paris
%        with influence factor greater than 5% of the most influential pair
%        will be included
% disflag - logical indication whether or not to display progress
%           information

% Created: June 17, 2015      Wenyu Li
% Modified: June 28, 2015      Wenyu Li    (sensitivity added)
% Modified: July 13, 2015     Wenyu Li   (Gradient and hessian provided)


frac = b2bopt.ExtraLinFraction;
disflag = b2bopt.Display; 
vars = obj.Variables.Values;
% LB = [vars.LowerBound]';
% UB = [vars.UpperBound]';
n_variable = length(vars);
LB = -ones(n_variable,1);
UB = -LB;
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end
flag = quadratictest(obj);
allVarnames = obj.VarNames;
if flag
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','Hessian',...
      'user-supplied','HessFcn',@hessianfcn,'TolFun',1e-10,'TolCon',1e-10);
else
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','MaxIter',1500,'TolFun',1e-10,'TolCon',1e-10);
end
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
n_inner = 1;
if n_inner > 1
   yin_result = zeros(n_inner,1);
   x0 = cell(n_inner,1);
   if disflag
      h = waitbar(0,'Start to calculate');
   end
   for j = 1:n_inner
      if isempty(obj.FeasiblePoint)
         [x0{j},yin_result(j)] = fmincon(@funxmin,[1;LB+(UB-LB).*lhsdesign(length(allVarnames),1)],...
            [],[],[],[],[-Inf;LB],[Inf;UB],@neq,opt);
      else
         [x0{j},yin_result(j)] = fmincon(@funxmin,[0;obj.FeasiblePoint],...
            [],[],[],[],[-Inf;LB],[Inf;UB],@neq,opt);
      end
      if disflag
         waitbar(j/n_inner,h,['Searching inner bound progress: ' int2str(j) ' / ' int2str(n_inner)])
      end
      yin_result(j) = -yin_result(j);
   end
   if disflag
      delete(h);
   end
else
   if isempty(obj.FeasiblePoint)
      [x0{1},yin_result] = fmincon(@funxmin,[1;LB+(UB-LB).*lhsdesign(length(allVarnames),1)],...
         [],[],[],[],[-Inf;LB],[Inf;UB],@neq,opt);
   else
      [x0{1},yin_result] = fmincon(@funxmin,[0;obj.FeasiblePoint],...
         [],[],[],[],[-Inf;LB],[Inf;UB],@neq,opt);
   end
   yin_result = -yin_result;
end
[yin,id] = max(yin_result);
if yin >= 0
   mx = 0.5*([vars.LowerBound]'+[vars.UpperBound]');
   dx = 0.5*([vars.UpperBound]'-[vars.LowerBound]');
   T = [1, zeros(1,n_variable);
      mx, diag(dx)];
   xtmp = T*x0{id};
   obj.FeasiblePoint = xtmp(2:end);
end

if disflag
   disp('=======================================================');
   disp('Calculating outer bound...');
   disp('=======================================================');
end
if flag
    [yout,sensitivity] = obj.sedumiconsisquadrel(b2bopt, abE);
%   [yout,sensitivity] = obj.cvxconsisquadrel(frac);
   obj.ConsistencyMeasure = [yin yout];
   obj.ConsistencySensitivity = sensitivity;
else
   n_outer = 1;
   yout_result = zeros(n_outer,1);
   sensitivity = cell(n_outer,1);
   if n_outer > 1
      if disflag
         h = waitbar(0,'Start to calculate');
      end
      for j = 1:n_outer
         [yout_result(j),sensitivity{j}] = obj.sedumiconsisrel(yin,frac,tol);
         if disflag
            waitbar(j/n_outer,h,['Progress: ' int2str(j) ' / ' int2str(n_outer)])
         end
      end
      if disflag
         delete(h)
      end
   else
      [yout_result(1),sensitivity{1}] = obj.sedumiconsisrel(yin,b2bopt,abE);
%      [yout_result(1),sensitivity{1}] = obj.cvxconsisrel(yin,frac,tol);
   end
   [yout,id] = min(yout_result);
   obj.ConsistencyMeasure = [yin yout];
   obj.ConsistencySensitivity = sensitivity{id};
end
if disflag
   disp(' ')
   disp('The calculation is done')
   disp(['Consistency LB: ' num2str(yin)])
   disp(['Consistency UB: ' num2str(yout)])
end



   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1;zeros(n_variable,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+1,2*n_units);
      for i = 1:n_units
         my = units(i).SurrogateModel.yScale.my;
         dy = units(i).SurrogateModel.yScale.dy;
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         l = (l-my)/dy;
         d = (d-my)/dy;
         u = (u-my)/dy;
         model = units(i).SurrogateModel;
         [~,id1,id2] = intersect(allVarnames,model.VarNames);
         id1 = id1+1;
         id3 = [1;id2+1];
         if isa(model,'B2BDC.B2Bmodels.RQModel')
            N = zeros(length(id3),length(id3));
            D = zeros(length(id3),length(id3));
            N1 = model.Numerator;
            D1 = model.Denominator;
            for i1 = 1:length(id3)
               for i2 = 1:length(id3)
                  N(i1,i2) = N1(id3(i1),id3(i2));
                  D(i1,i2) = D1(id3(i1),id3(i2));
               end
            end
            c(2*i-1,1) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)])-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(1) = u-d;
            grad1(id1) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad2(1) = d-l;
            grad2(id1) = -2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(model,'B2BDC.B2Bmodels.QModel')
            quadCoef = zeros(length(id3),length(id3));
            quadCoef2 = model.CoefMatrix;
            for i1 = 1:length(id3)
               for i2 = 1:length(id3)
                  quadCoef(i1,i2) = quadCoef2(id3(i1),id3(i2));
               end
            end
            c(2*i-1,1) = [1;x(id1)]'*quadCoef*[1;x(id1)]-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-[1;x(id1)]'*quadCoef*[1;x(id1)];
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            grad1(1) = u-d;
            grad2(id1) = -2*quadCoef(2:end,2:end)*x(id1)-2*quadCoef(2:end,1);
            grad2(1) = d-l;
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      ceq = [];
      geq = [];
   end

   function hessian = hessianfcn(x, lambda)
      hessian = zeros(n_variable+1,n_variable+1);
      for i = 1:n_units
         model = units(i).SurrogateModel;
         modelVar = model.VarNames;
         [~,id1,id2] = intersect(allVarnames,modelVar);
         id1 = id1+1;
         quadHessian = model.Hessian;
         hessianMatrix1 = zeros(n_variable+1,n_variable+1);
         hessianMatrix2 = zeros(n_variable+1,n_variable+1);
         for i2 = 1:length(id2)
            for i3 = 1:length(id2)
               hessianMatrix1(id1(i2),id1(i3)) = quadHessian(id2(i2),id2(i3));
               hessianMatrix2(id1(i2),id1(i3)) = -quadHessian(id2(i2),id2(i3));
            end
         end
         hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
            lambda.ineqnonlin(2*i)*hessianMatrix2;
      end
   end

end