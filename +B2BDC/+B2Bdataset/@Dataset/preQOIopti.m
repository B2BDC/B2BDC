<<<<<<< Updated upstream
function [Qmin, Qmax, s, xOpt, abE] = preQOIopti(obj,QOIobj,disflag,rflag,b2bopt)
% subfunction to calculate QOI inner bound with opti

%  Created: Oct 10, 2016     Wenyu Li

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
   n1 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub-Alb)./sum(abs(A0),2);
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
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
allVarnames = obj.VarNames;
model0 = QOIobj;
[~,~,id0] = intersect(model0.VarNames,allVarnames,'stable');
if isa(model0,'B2BDC.B2Bmodels.RQModel')
   N0 = model0.Numerator;
   D0 = model0.Denominator;
elseif isa(model0,'B2BDC.B2Bmodels.QModel')
   Coef0 = model0.CoefMatrix;
end
units = obj.DatasetUnits.Values;
if rflag
   n_units = length(units)-1;
else
   n_units = length(units);
end
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
x0 = obj.FeasiblePoint;
Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,Qmin,exitflag,info] = solve(Opt);
xOpt(:,1) = x;
lam = info.Lambda;
s.min.expu = -lam.ineqnonlin(1:2:2*n_units).*bds;
s.min.expl = -lam.ineqnonlin(2:2:2*n_units).*bds;
s.min.varu = -lam.upper.*xbd;
s.min.varl = -lam.lower.*xbd;
if ~isempty(A1)
   s.min.linu = zeros(n1,1);
   s.min.linl = zeros(n1,1);
   s.min.linu(~ieq) = -lam.ineqlin(1:n1).*dA(~ieq);
   s.min.linl(~ieq) = -lam.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.min.linu = [];
   s.min.linl = [];
end
Opt = opti('fun',@funxmax,'grad',@gradmax,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,Qmax,exitflag,info] = solve(Opt);
Qmax = -Qmax;
xOpt(:,2) = x;
lam = info.Lambda;
s.max.expu = lam.ineqnonlin(1:2:2*n_units).*bds;
s.max.expl = lam.ineqnonlin(2:2:2*n_units).*bds;
s.max.varu = lam.upper.*xbd;
s.max.varl = lam.lower.*xbd;
if ~isempty(A1)
   s.max.linu = zeros(n1,1);
   s.max.linl = zeros(n1,1);
   s.max.linu(~ieq) = lam.ineqlin(1:n1).*dA(~ieq);
   s.max.linl(~ieq) = lam.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.max.linu = [];
   s.max.linl = [];
end



   function y = funxmin(x)
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = ([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = [1;x(id0)]'*Coef0*[1;x(id0)];
      end
   end

   function gy = gradmin(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         gy(id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         gy(id0) = 2*Coef0(2:end,2:end)*x(id0)+2*Coef0(2:end,1);
      end
      gy = gy';
   end

   function y = funxmax(x)
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = -([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = -[1;x(id0)]'*Coef0*[1;x(id0)];
      end
   end

   function gy = gradmax(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         gy(id0) = -2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         gy(id0) = -2*Coef0(2:end,2:end)*x(id0)-2*Coef0(2:end,1);
      end
      gy = gy';
   end

   function c = nlcon(x)
      c = zeros(2*n_units,1);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)])-u;
            c(2*i,1) = l-([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id)]'*quadCoef*[1;x(id)]-u;
            c(2*i,1) =  l-[1;x(id)]'*quadCoef*[1;x(id)];
         end
      end
   end

   function g = jac(x)
      g = zeros(n_variable,2*n_units);
      for i = 1:n_units
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad2(id) = -2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*quadCoef(2:end,2:end)*x(id)+2*quadCoef(2:end,1);
            grad2(id) = -2*quadCoef(2:end,2:end)*x(id)-2*quadCoef(2:end,1);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      g = g';
   end

   function hessian = H(x, sigma,lambda)
      h0 = zeros(n_variable);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         h0(id0,id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.RQModel')
         h0(id0,id0) = model0.Hessian;
      end
      hessian = sigma*h0;
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
end
=======
function [Qmin, Qmax, s, xOpt, abE] = preQOIopti(obj,QOIobj,disflag,rflag,b2bopt)
% subfunction to calculate QOI inner bound with opti

%  Created: Oct 10, 2016     Wenyu Li

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
   n1 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub-Alb)./sum(abs(A0),2);
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
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
allVarnames = obj.VarNames;
model0 = QOIobj;
[~,~,id0] = intersect(model0.VarNames,allVarnames,'stable');
if isa(model0,'B2BDC.B2Bmodels.RQModel')
   N0 = model0.Numerator;
   D0 = model0.Denominator;
elseif isa(model0,'B2BDC.B2Bmodels.QModel')
   Coef0 = model0.CoefMatrix;
end
units = obj.DatasetUnits.Values;
if rflag
   n_units = length(units)-1;
else
   n_units = length(units);
end
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
x0 = obj.FeasiblePoint;
Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,Qmin,exitflag,info] = solve(Opt);
xOpt(:,1) = x;
lam = info.Lambda;
s.min.expu = -lam.ineqnonlin(1:2:2*n_units).*bds;
s.min.expl = -lam.ineqnonlin(2:2:2*n_units).*bds;
s.min.varu = -lam.upper.*xbd;
s.min.varl = -lam.lower.*xbd;
if ~isempty(A1)
   s.min.linu = zeros(n1,1);
   s.min.linl = zeros(n1,1);
   s.min.linu(~ieq) = -lam.ineqlin(1:n1).*dA(~ieq);
   s.min.linl(~ieq) = -lam.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.min.linu = [];
   s.min.linl = [];
end
Opt = opti('fun',@funxmax,'grad',@gradmax,'nlmix',@nlcon,zeros(2*n_units,1),...
   -ones(2*n_units,1),'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[x,Qmax,exitflag,info] = solve(Opt);
Qmax = -Qmax;
xOpt(:,2) = x;
lam = info.Lambda;
s.max.expu = lam.ineqnonlin(1:2:2*n_units).*bds;
s.max.expl = lam.ineqnonlin(2:2:2*n_units).*bds;
s.max.varu = lam.upper.*xbd;
s.max.varl = lam.lower.*xbd;
if ~isempty(A1)
   s.max.linu = zeros(n1,1);
   s.max.linl = zeros(n1,1);
   s.max.linu(~ieq) = lam.ineqlin(1:n1).*dA(~ieq);
   s.max.linl(~ieq) = lam.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.max.linu = [];
   s.max.linl = [];
end



   function y = funxmin(x)
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = ([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = [1;x(id0)]'*Coef0*[1;x(id0)];
      end
   end

   function gy = gradmin(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         gy(id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         gy(id0) = 2*Coef0(2:end,2:end)*x(id0)+2*Coef0(2:end,1);
      end
      gy = gy';
   end

   function y = funxmax(x)
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = -([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = -[1;x(id0)]'*Coef0*[1;x(id0)];
      end
   end

   function gy = gradmax(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         gy(id0) = -2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         gy(id0) = -2*Coef0(2:end,2:end)*x(id0)-2*Coef0(2:end,1);
      end
      gy = gy';
   end

   function c = nlcon(x)
      c = zeros(2*n_units,1);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)])-u;
            c(2*i,1) = l-([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id)]'*quadCoef*[1;x(id)]-u;
            c(2*i,1) =  l-[1;x(id)]'*quadCoef*[1;x(id)];
         end
      end
   end

   function g = jac(x)
      g = zeros(n_variable,2*n_units);
      for i = 1:n_units
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad2(id) = -2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*quadCoef(2:end,2:end)*x(id)+2*quadCoef(2:end,1);
            grad2(id) = -2*quadCoef(2:end,2:end)*x(id)-2*quadCoef(2:end,1);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      g = g';
   end

   function hessian = H(x, sigma,lambda)
      h0 = zeros(n_variable);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         h0(id0,id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.RQModel')
         h0(id0,id0) = model0.Hessian;
      end
      hessian = sigma*h0;
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
end
>>>>>>> Stashed changes
