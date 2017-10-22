<<<<<<< Updated upstream
function [yin_result,s,xopt,abE,flag] = relCMopti(obj,disflag,b2bopt)
% subfunction to calculate CM inner bound with fmincon

%  Created: Oct 7, 2016     Wenyu Li

if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
xbd = UB-LB;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   n3 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub - Alb)./sum(abs(A0),2);
   dA = Aub - Alb;
   ieq = (eqTest <= tolerance);
   if any(ieq)
      A1 = [A0(~ieq,:); -A0(~ieq,:)];
      B1 = [Aub(~ieq); -Alb(~ieq)];
      Aeq = A0(ieq,:);
      Beq = 0.5*(Alb(ieq)+Aub(ieq));
   else
      A1 = [A0; -A0];
      Aeq = [];
      B1 = [Aub; -Alb];
      Beq = [];
   end
   n1 = size(A1,1);
   n2 = size(Aeq,1);
   A1 = [zeros(n1,1), A1];
   if n2 ~= 0
      Aeq = [zeros(n2,1), Aeq];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
bds = obj.calBound;
bds = diff(bds,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end
flag = quadratictest(obj);
allVarnames = obj.VarNames;
idall = cell(n_units,1);
Nall = cell(n_units,1);
Dall = cell(n_units,1);
Qall = cell(n_units,1);
for j = 1:n_units
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
      Nall{j} = tmodel.Numerator;
      Dall{j} = tmodel.Denominator;
   elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
      Qall{j} = tmodel.CoefMatrix;
   end
end
optSet = optiset;
optSet.tolafun = 1e-6;
if isempty(obj.FeasiblePoint)
   x0 = [0;LB+(UB-LB).*lhsdesign(n_variable,1)];
else
   x0 = [0;obj.FeasiblePoint];
end
% Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
%    -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
%    'jac',@jac,'ndec',2,'hess',@H);
Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,fval,exitflag,info] = solve(Opt);
lam = info.Lambda;
s.expu = lam.ineqnonlin(1:2:2*n_units).*bds;
s.expl = lam.ineqnonlin(2:2:2*n_units).*bds;
s.varu = lam.upper(2:end).*xbd;
s.varl = lam.lower(2:end).*xbd;
if ~isempty(A1)
   s.linu = zeros(n3,1);
   s.linl = zeros(n3,1);
   s.linu(~ieq) = lam.ineqlin(1:n3).*dA(~ieq);
   s.linl(~ieq) = lam.ineqlin(n3+1:2*n3).*dA(~ieq);
else
   s.linu = [];
   s.linl = [];
end
yin_result = -fval;
xopt = x(2:end);

   function y = funxmin(x)
      y = -x(1);
   end

   function gy = gradmin(x)
      gy = [-1;zeros(n_variable,1)]';
   end

   function c = nlcon(x)
      c = zeros(2*n_units,1);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)])-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id+1)]'*quadCoef*[1;x(id+1)]-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-[1;x(id+1)]'*quadCoef*[1;x(id+1)];
         end
      end
   end

   function g = jac(x)
      g = zeros(n_variable+1,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(1) = u-d;
            grad1(id+1) = 2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            grad2(1) = d-l;
            grad2(id+1) = -2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(id+1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
            grad1(1) = u-d;
            grad2(id+1) = -2*quadCoef(2:end,2:end)*x(id+1)-2*quadCoef(2:end,1);
            grad2(1) = d-l;
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      g = g';
   end

   function hessian = H(x, sigma,lambda)
      hessian = sigma*zeros(n_variable+1);
      for i = 1:n_units
         id = idall{i};
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable+1);
         hessianMatrix2 = zeros(n_variable+1);
         hessianMatrix1(id+1,id+1) = quadHessian;
         hessianMatrix2(id+1,id+1) = -quadHessian;
         hessian = hessian+lambda(2*i-1)*(tril(hessianMatrix1)+tril(hessianMatrix1,-1))+...
            lambda(2*i)*(tril(hessianMatrix2)+tril(hessianMatrix2,-1));
      end
      hessian = sparse(hessian);
   end
=======
function [yin_result,s,xopt,abE,flag] = relCMopti(obj,disflag,b2bopt)
% subfunction to calculate CM inner bound with fmincon

%  Created: Oct 7, 2016     Wenyu Li

if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
xbd = UB-LB;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   n3 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub - Alb)./sum(abs(A0),2);
   dA = Aub - Alb;
   ieq = (eqTest <= tolerance);
   if any(ieq)
      A1 = [A0(~ieq,:); -A0(~ieq,:)];
      B1 = [Aub(~ieq); -Alb(~ieq)];
      Aeq = A0(ieq,:);
      Beq = 0.5*(Alb(ieq)+Aub(ieq));
   else
      A1 = [A0; -A0];
      Aeq = [];
      B1 = [Aub; -Alb];
      Beq = [];
   end
   n1 = size(A1,1);
   n2 = size(Aeq,1);
   A1 = [zeros(n1,1), A1];
   if n2 ~= 0
      Aeq = [zeros(n2,1), Aeq];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
bds = obj.calBound;
bds = diff(bds,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end
flag = quadratictest(obj);
allVarnames = obj.VarNames;
idall = cell(n_units,1);
Nall = cell(n_units,1);
Dall = cell(n_units,1);
Qall = cell(n_units,1);
for j = 1:n_units
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
      Nall{j} = tmodel.Numerator;
      Dall{j} = tmodel.Denominator;
   elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
      Qall{j} = tmodel.CoefMatrix;
   end
end
optSet = optiset;
optSet.tolafun = 1e-6;
if isempty(obj.FeasiblePoint)
   x0 = [0;LB+(UB-LB).*lhsdesign(n_variable,1)];
else
   x0 = [0;obj.FeasiblePoint];
end
% Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
%    -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
%    'jac',@jac,'ndec',2,'hess',@H);
Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,fval,exitflag,info] = solve(Opt);
lam = info.Lambda;
s.expu = lam.ineqnonlin(1:2:2*n_units).*bds;
s.expl = lam.ineqnonlin(2:2:2*n_units).*bds;
s.varu = lam.upper(2:end).*xbd;
s.varl = lam.lower(2:end).*xbd;
if ~isempty(A1)
   s.linu = zeros(n3,1);
   s.linl = zeros(n3,1);
   s.linu(~ieq) = lam.ineqlin(1:n3).*dA(~ieq);
   s.linl(~ieq) = lam.ineqlin(n3+1:2*n3).*dA(~ieq);
else
   s.linu = [];
   s.linl = [];
end
yin_result = -fval;
xopt = x(2:end);

   function y = funxmin(x)
      y = -x(1);
   end

   function gy = gradmin(x)
      gy = [-1;zeros(n_variable,1)]';
   end

   function c = nlcon(x)
      c = zeros(2*n_units,1);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)])-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id+1)]'*quadCoef*[1;x(id+1)]-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-[1;x(id+1)]'*quadCoef*[1;x(id+1)];
         end
      end
   end

   function g = jac(x)
      g = zeros(n_variable+1,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(1) = u-d;
            grad1(id+1) = 2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            grad2(1) = d-l;
            grad2(id+1) = -2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(id+1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
            grad1(1) = u-d;
            grad2(id+1) = -2*quadCoef(2:end,2:end)*x(id+1)-2*quadCoef(2:end,1);
            grad2(1) = d-l;
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      g = g';
   end

   function hessian = H(x, sigma,lambda)
      hessian = sigma*zeros(n_variable+1);
      for i = 1:n_units
         id = idall{i};
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable+1);
         hessianMatrix2 = zeros(n_variable+1);
         hessianMatrix1(id+1,id+1) = quadHessian;
         hessianMatrix2(id+1,id+1) = -quadHessian;
         hessian = hessian+lambda(2*i-1)*(tril(hessianMatrix1)+tril(hessianMatrix1,-1))+...
            lambda(2*i)*(tril(hessianMatrix2)+tril(hessianMatrix2,-1));
      end
      hessian = sparse(hessian);
   end
>>>>>>> Stashed changes
end