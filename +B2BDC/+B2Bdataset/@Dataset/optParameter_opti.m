<<<<<<< Updated upstream
function xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,w,b2bopt)

% function XOPT = OPTPARAMETER_OPTI(OBJ,ID,FLAG,NFLAG,HFLAG,W,B2BOPT) 
% returns the optimal parameter according to the optimization standard

%  Created: Nov 15, 2016    Wenyu Li

vList = obj.Variables;
nVar = vList.Length;
nomVal = [vList.Values.NominalValue]';
xbd = vList.calBound;
LB = xbd(:,1);
UB = xbd(:,2);
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub - Alb)./sum(abs(A0),2);
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
if isempty(obj.FeasiblePoint)
   x0 = vList.makeLHSsample(1)';
else
   x0 = obj.FeasiblePoint;
end
if ~isempty(id)
   nt = length(id);
   tmpA = zeros(nt,nVar);
   tmpA(:,id) = eye(nt);
   tmpB = nomVal(id);
   Aeq = [Aeq;tmpA];
   Beq = [Beq;tmpB];
   x0(id) = 0;
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
qbd = obj.calBound;
bds = diff(qbd,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
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
% Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
%    -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
%    'jac',@jac,'ndec',2,'hess',@H);
Opt = opti('fun',@funxmin,'grad',@gradmin,'nl',@nlcon,qbd(:,1)-abE,...
   qbd(:,2)+abE,'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[xopt,yopt,exitflag,info] = solve(Opt);
if exitflag < 0
   xopt = [];
end






   function y = funxmin(x)
      y = 0;
      if ~Nflag
         for i = 1:n_units
            tw = w(i);
            id1 = idall{i};
            d = units(i).ObservedValue;
            if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
               N = Nall{i};
               D = Dall{i};
               tmpN = [1;x(id1)]'*N*[1;x(id1)];
               tmpD = [1;x(id1)]'*D*[1;x(id1)];
               y = y+tw^2*(tmpN/tmpD-d)^2;
            elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
               quadCoef = Qall{i};
               tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
               y = y + tw^2*(tmpQ - d)^2;
            end
         end
      else
         y = norm(x.*w,1);
      end
   end

   function gy = gradmin(x)
      gy = zeros(nVar,1);
      if ~Nflag
         for i = 1:n_units
            tw = w(i);
            id1 = idall{i};
            d = units(i).ObservedValue;
            if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
               N = Nall{i};
               D = Dall{i};
               tmpN = [1;x(id1)]'*N*[1;x(id1)];
               tmpD = [1;x(id1)]'*D*[1;x(id1)];
               gy(id1) = gy(id1)+2*tw^2*(tmpN/tmpD-d)/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
               quadCoef = Qall{i};
               tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
               gy(id1) = gy(id1) + 2*tw^2*(tmpQ-d)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      else
         gy(x~=0) = abs(w).*sign(x(x~=0));
      end
      gy = gy';
   end

   function c = nlcon(x)
      c = zeros(n_units,1);
      for i = 1:n_units
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(i,1) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(i,1) = [1;x(id1)]'*quadCoef*[1;x(id1)];
         end
      end
   end
         
   function g = jac(x)
      g = zeros(n_units,nVar);
      for i = 1:n_units
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad = zeros(1,nVar);
            grad(id1) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
         end
         g(i,:) = grad;
      end
   end

   function hessian = H(x, sigma,lambda)
      hessian = sigma*zeros(nVar);
      for i = 1:n_units
         id1 = idall{i};
         tw = w(i);
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable);
         hessianMatrix1(id1,id1) = quadHessian;
         hessianMatrix2 = 2*tw^2*hessianMatrix1;
         hessian = hessian+lambda(i)*(tril(hessianMatrix1)+tril(hessianMatrix1,-1))+...
            sigma*(tril(hessianMatrix2)+tril(hessianMatrix2,-1));
      end
      hessian = sparse(hessian);
   end
         
         
         
         
         
         
=======
function xopt = optParameter_opti(obj,id,flag,Nflag,Hflag,w,b2bopt)

% function XOPT = OPTPARAMETER_OPTI(OBJ,ID,FLAG,NFLAG,HFLAG,W,B2BOPT) 
% returns the optimal parameter according to the optimization standard

%  Created: Nov 15, 2016    Wenyu Li

vList = obj.Variables;
nVar = vList.Length;
nomVal = [vList.Values.NominalValue]';
xbd = vList.calBound;
LB = xbd(:,1);
UB = xbd(:,2);
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub - Alb)./sum(abs(A0),2);
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
if isempty(obj.FeasiblePoint)
   x0 = vList.makeLHSsample(1)';
else
   x0 = obj.FeasiblePoint;
end
if ~isempty(id)
   nt = length(id);
   tmpA = zeros(nt,nVar);
   tmpA(:,id) = eye(nt);
   tmpB = nomVal(id);
   Aeq = [Aeq;tmpA];
   Beq = [Beq;tmpB];
   x0(id) = 0;
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
qbd = obj.calBound;
bds = diff(qbd,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
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
% Opt = opti('fun',@funxmin,'grad',@gradmin,'nlmix',@nlcon,zeros(2*n_units,1),...
%    -ones(2*n_units,1),'bounds',[-Inf;LB],[1;UB],'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
%    'jac',@jac,'ndec',2,'hess',@H);
Opt = opti('fun',@funxmin,'grad',@gradmin,'nl',@nlcon,qbd(:,1)-abE,...
   qbd(:,2)+abE,'bounds',LB,UB,'x0',x0,'ineq',A1,B1,'eq',Aeq,Beq,...
   'jac',@jac,'ndec',2,'options',optSet);
[xopt,yopt,exitflag,info] = solve(Opt);
if exitflag < 0
   xopt = [];
end






   function y = funxmin(x)
      y = 0;
      if ~Nflag
         for i = 1:n_units
            tw = w(i);
            id1 = idall{i};
            d = units(i).ObservedValue;
            if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
               N = Nall{i};
               D = Dall{i};
               tmpN = [1;x(id1)]'*N*[1;x(id1)];
               tmpD = [1;x(id1)]'*D*[1;x(id1)];
               y = y+tw^2*(tmpN/tmpD-d)^2;
            elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
               quadCoef = Qall{i};
               tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
               y = y + tw^2*(tmpQ - d)^2;
            end
         end
      else
         y = norm(x.*w,1);
      end
   end

   function gy = gradmin(x)
      gy = zeros(nVar,1);
      if ~Nflag
         for i = 1:n_units
            tw = w(i);
            id1 = idall{i};
            d = units(i).ObservedValue;
            if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
               N = Nall{i};
               D = Dall{i};
               tmpN = [1;x(id1)]'*N*[1;x(id1)];
               tmpD = [1;x(id1)]'*D*[1;x(id1)];
               gy(id1) = gy(id1)+2*tw^2*(tmpN/tmpD-d)/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
               quadCoef = Qall{i};
               tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
               gy(id1) = gy(id1) + 2*tw^2*(tmpQ-d)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      else
         gy(x~=0) = abs(w).*sign(x(x~=0));
      end
      gy = gy';
   end

   function c = nlcon(x)
      c = zeros(n_units,1);
      for i = 1:n_units
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(i,1) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(i,1) = [1;x(id1)]'*quadCoef*[1;x(id1)];
         end
      end
   end
         
   function g = jac(x)
      g = zeros(n_units,nVar);
      for i = 1:n_units
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            grad = zeros(1,nVar);
            grad(id1) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            grad(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
         end
         g(i,:) = grad;
      end
   end

   function hessian = H(x, sigma,lambda)
      hessian = sigma*zeros(nVar);
      for i = 1:n_units
         id1 = idall{i};
         tw = w(i);
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable);
         hessianMatrix1(id1,id1) = quadHessian;
         hessianMatrix2 = 2*tw^2*hessianMatrix1;
         hessian = hessian+lambda(i)*(tril(hessianMatrix1)+tril(hessianMatrix1,-1))+...
            sigma*(tril(hessianMatrix2)+tril(hessianMatrix2,-1));
      end
      hessian = sparse(hessian);
   end
         
         
         
         
         
         
>>>>>>> Stashed changes
end