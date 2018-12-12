function [d,xopt] = calculateDistanceCurve(obj,x0,b2bopt,C)
% D=CALCULATEDISTANCECURVE(OBJ,X,OPT,C) calculates the minimum average difference
% between the speicified cureve and any curve in the feasible set.

%  Created: Sep 17, 2018     Wenyu Li

if size(x0,2) == 1
   x0 = x0';
end
starget = C.Value;
funTrue = C.function;
ntarget = length(starget);
xtarget = funTrue(starget)+repmat(x0,ntarget,1);
vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
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
xFea = obj.FeasiblePoint;
if obj.ModelDiscrepancyFlag
   MDvar = obj.ModelDiscrepancy.Variables;
   nMD = MDvar.Length;
   MDbd = MDvar.calBound;
   LB = [LB; MDbd(:,1)];
   UB = [UB; MDbd(:,2)];
else
   nMD = 0;
end
if obj.ParameterDiscrepancyFlag
   PDvar = obj.ParameterDiscrepancy.Variables;
   nPD = PDvar.Length;
   PDbd = PDvar.calBound;
   npd = obj.ParameterDiscrepancy.CorrectionDimension;
   pdBasis = cell(n_variable,1);
   fPD = obj.ParameterDiscrepancy.BasisFunction;
   for j = 1:n_variable
      if npd(j) ~= 0
         pdBasis{j} = fPD{1}(starget);
         fPD(1) = [];
      end
   end
   LB = [LB; PDbd(:,1)];
   UB = [UB; PDbd(:,2)];
else
   npd = zeros(1,n_variable);
   nPD = 0;
end
if nMD > 0
   xFea = [xFea; obj.ModelDiscrepancy.FeasiblePoint];
end
if nPD > 0
   xFea = [xFea; obj.ParameterDiscrepancy.FeasiblePoint];
end
[idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion;
if ~isempty(A1)
   A1 = [A1 zeros(size(A1,1),nMD+nPD)]; 
end
if ~isempty(Aeq)
   Aeq = [Aeq zeros(size(Aeq,1),nMD+nPD)];
end
if ~isempty(APD)
   infFlag = bPD == inf;
   APD(infFlag,:) = [];
   bPD(infFlag) = [];
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','MaxIter',3000,...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
nStart = b2bopt.LocalStart;
b2bopt.SampleOption.StepInterval = 100;
xSample = obj.collectHRsamples(nStart-1);
xStart = [xFea'; xSample.x];
xmin = cell(nStart,1);
ymin = zeros(nStart,1);
minflag = zeros(nStart,1);

for j = 1:nStart
   [xmin{j},ymin(j),minflag(j)] = fmincon(@funxmin,xStart(j,:)',...
      [A1;APD],[B1;bPD],Aeq,Beq,LB,UB,@neq,opt);
end
xmin(minflag<=0) = [];
ymin(minflag<=0) = [];
if ~isempty(ymin)
   [dd,id] = min(ymin);
   d = sqrt(dd);
   xopt.x = xmin{id};
   xopt.dimension = [n_variable nMD nPD];
else
   d = inf;
   xopt = [];
end


   function [y,gy] = funxmin(x)
      gy = zeros(n_variable+nMD+nPD,1);
      dx = zeros(ntarget,n_variable);
      for i = 1:n_variable
         if npd(i) ~= 0
            dx(:,i) = pdBasis{i}*x(n_variable+nMD+sum(npd(1:i-1))+1:n_variable+nMD+sum(npd(1:i)));
         end
      end
      xx = repmat(x(1:n_variable)',ntarget,1);
      y = xx+dx-xtarget;
      gy(1:n_variable) = 2*mean(y);
      for i = 1:n_variable
         if npd(i) ~= 0
            gy(n_variable+nMD+sum(npd(1:i-1))+1:n_variable+nMD+sum(npd(1:i))) = ...
               mean(repmat(2*y(:,i),1,npd(i)).*pdBasis{i});
         end
      end
      y = mean(sum(y.^2,2));
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD,2*n_units);
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
            g(id1,2*i-1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            g(id1,2*i) = -g(id1,2*i-1);
         end
      end
      ceq = [];
      geq = [];
   end
end