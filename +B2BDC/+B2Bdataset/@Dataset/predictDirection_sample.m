function uq = predictDirection_sample(obj,vv,nPC,xAve,d)

vPC = vv(:,1:nPC);
veq = vv(:,nPC+1:end);
y0 = xAve*vv;
dd = vPC*d;
f0 = xAve*dd;
vars = obj.Variables.Values;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
allVarnames = obj.VarNames;
model0 = generateModel([0;dd],vars);
[~,~,id0] = intersect(model0.VarNames,allVarnames,'stable');
if isa(model0,'B2BDC.B2Bmodels.RQModel')
   N0 = model0.Numerator;
   D0 = model0.Denominator;
elseif isa(model0,'B2BDC.B2Bmodels.QModel')
   Coef0 = model0.CoefMatrix;
end
units = obj.DatasetUnits.Values;
abE = zeros(length(units),1);
bds = obj.calBound;
bds = diff(bds,[],2);
n_units = length(units);
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
[~,Qmin,ttt,~,ilam_min] = fmincon(@funxmin,obj.FeasiblePoint,[],[],veq',...
   y0(nPC+1:end)',LB,UB,@neq,opt);
[~,Qmax,ttt,~,ilam_max] = fmincon(@funxmax,obj.FeasiblePoint,[],[],veq',...
   y0(nPC+1:end)',LB,UB,@neq,opt);
Qmax = -Qmax;
uq = [Qmin-f0, Qmax-f0];






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