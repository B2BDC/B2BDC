function [Qmin, Qmax, s, xOpt, abE] = preQOIfmincon(obj,QOIobj,disflag,rflag,b2bopt,xs_min,xs_max)
% subfunction to calculate QOI inner bound with fmincon

%  Created: Oct 10, 2016     Wenyu Li

vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
dx = UB-LB;
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
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point',...
   'TolFun',1e-6,'TolCon',1e-6);
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
[xmin,Qmin,ttt,~,ilam_min] = fmincon(@funxmin,xs_min,A1,B1,Aeq,...
   Beq,LB,UB,@neq,opt);
[xmax,Qmax,ttt,~,ilam_max] = fmincon(@funxmax,xs_max,A1,B1,Aeq,...
   Beq,LB,UB,@neq,opt);
Qmax = -Qmax;
s.min.expu = -ilam_min.ineqnonlin(1:2:2*n_units).*bds;
s.min.expl = -ilam_min.ineqnonlin(2:2:2*n_units).*bds;
s.min.varu = -ilam_min.upper.*dx;
s.min.varl = -ilam_min.lower.*dx;
if ~isempty(A1)
   s.min.linu = zeros(n1,1);
   s.min.linl = zeros(n1,1);
   s.min.linu(~ieq) = -ilam_min.ineqlin(1:n1).*dA(~ieq);
   s.min.linl(~ieq) = -ilam_min.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.min.linu = [];
   s.min.linl = [];
end

s.max.expu = ilam_max.ineqnonlin(1:2:2*n_units).*bds;
s.max.expl = ilam_max.ineqnonlin(2:2:2*n_units).*bds;
s.max.varu = ilam_max.upper.*dx;
s.max.varl = ilam_max.lower.*dx;
if ~isempty(A1)
   s.max.linu = zeros(n1,1);
   s.max.linl = zeros(n1,1);
   s.max.linu(~ieq) = ilam_max.ineqlin(1:n1).*dA(~ieq);
   s.max.linl(~ieq) = ilam_max.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.max.linu = [];
   s.max.linl = [];
end

xOpt(:,1) = xmin;
xOpt(:,2) = xmax;





   function [y,gy] = funxmin(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = ([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
         gy(id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = [1;x(id0)]'*Coef0*[1;x(id0)];
         gy(id0) = 2*Coef0(2:end,2:end)*x(id0)+2*Coef0(2:end,1);
      end
   end

   function [y,gy] = funxmax(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = -([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
         gy(id0) = -2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = -[1;x(id0)]'*Coef0*[1;x(id0)];
         gy(id0) = -2*Coef0(2:end,2:end)*x(id0)-2*Coef0(2:end,1);
      end
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)])-u;
            c(2*i,1) =  l-([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id1) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad2(id1) = -2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id1)]'*quadCoef*[1;x(id1)]-u;
            c(2*i,1) =  l-[1;x(id1)]'*quadCoef*[1;x(id1)];
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            grad2(id1) = -2*quadCoef(2:end,2:end)*x(id1)-2*quadCoef(2:end,1);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      ceq = [];
      geq = [];
   end

   function hessian = hessianfcn(x, lambda)
      hessian = zeros(n_variable);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         hessian(id0,id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.RQModel')
         hessian(id0,id0) = model0.Hessian;
      end
      for i = 1:n_units
         id = idall{i};
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable);
         hessianMatrix2 = zeros(n_variable);
         hessianMatrix1(id,id) = quadHessian;
         hessianMatrix2(id,id) = -quadHessian;
         hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
            lambda.ineqnonlin(2*i)*hessianMatrix2;
      end
   end



end