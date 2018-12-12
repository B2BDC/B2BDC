function [xopt,yopt] = optimize_Parameter(obj,logFlag,opt,weight)
% XOPT=OPTIMIZE_PARAMETER(OBJ,LOGFLAG,OPT,WEIGHT) finds the optimal parameter
% subject to the optimization criteria given in opt. 

%  Created: Oct 25, 2017

if nargin < 2 || isempty(logFlag)
   logFlag = [true(obj.Length,1); false(obj.Variables.Length,1)];
end
if nargin < 3
   opt = generateOpt;
end
optimOpt = opt.OptimOption;
method = optimOpt.OptimizationMethod;
vList = obj.Variables;
ns = optimOpt.RandomStart;
if ~strcmp(method,'LSH')
   if ~obj.isConsistent(opt)
      disp('The dataset is inconsistent!')
      xopt = [];
      return
   else
      x0 = obj.collectSamples(ns,[],opt);
   end
   if strcmp(method,'1NF')
      optimOpt.PenaltyWeight = 'user-defined';
   end
else
%    x0 = vList.makeLHSsample(ns);
   x0 = obj.collectSamples(ns,[],opt);
end
nVar = vList.Length;
varNom = [vList.Values.NominalValue]';
H = vList.calBound;
LB = H(:,1);
UB = H(:,2);
fminopt = optimoptions('fmincon','Display','none','GradObj','on','MaxFunctionEvaluations',3000,...
      'GradConstr','on','Algorithm','interior-point','MaxIter',3000,'TolFun',1e-10,'TolCon',1e-10,...
      'StepTolerance',1e-20);
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
bds = obj.calBound;
obs = obj.calObserve;
units = obj.DatasetUnits.Values;
n_units = length(units);
obs(logFlag(1:n_units)) = exp(obs(logFlag(1:n_units)));
varNom(logFlag(n_units+1:end)) = exp(varNom(logFlag(n_units+1:end)));
abE = zeros(length(units),1);
if opt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
switch optimOpt.PenaltyWeight
   case 'relative'
      w = [1./obs; zeros(nVar,1)];
   case 'absolute'
      w = [ones(n_units,1); zeros(nVar,1)];
   case 'user-defined'
      if nargin > 3
         w = weight;
      elseif strcmp(optimOpt.OptimizationMethod,'1NF')
         w = [zeros(n_units,1); ones(nVar,1)];
      else
         w = [1./obs; zeros(nVar,1)];
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

xopt = zeros(nVar,ns);
yopt = zeros(1,ns);
exitflag = zeros(1,ns);
for j = 1:ns
   switch method
      case 'LSF'
         [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS,x0(j,:)',...
            A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
      case 'LSH'
         [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS,x0(j,:)',...
            A1,B1,Aeq,Beq,LB,UB,[],fminopt);
      case '1NF'
         [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_1N,x0(j,:)',...
            A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
   end
end
id = find(exitflag<0);
xopt(:,id) = [];
yopt(id) = [];
if isempty(yopt)
   disp('Optimal point is not found.')
else
   [~,id] = min(yopt);
   xopt = xopt(:,id);
   yopt = yopt(id);
end

   function [y,gy] = min_LS(x)
      y = 0;
      gy = zeros(nVar,1);
      unitID = find(w(1:n_units));
      for i = unitID'
         tw = w(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            if logFlag(i)
               y = y+tw^2*(exp(tmpN/tmpD)-obs(i))^2;
            else
               y = y+tw^2*(tmpN/tmpD-obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1)+2*tw^2*(exp(tmpN/tmpD)-obs(i))*exp(tmpN/tmpD)...
                  /tmpD^2*(tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            else
               gy(id1) = gy(id1)+2*tw^2*(tmpN/tmpD-obs(i))/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            end
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            if logFlag(i)
               y = y + tw^2*(exp(tmpQ) - obs(i))^2;
            else
               y = y + tw^2*(tmpQ - obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1) + 2*tw^2*(exp(tmpQ)-obs(i))*exp(tmpQ)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            else
               gy(id1) = gy(id1) + 2*tw^2*(tmpQ-obs(i))*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      end
      if strcmp(optimOpt.PenaltyWeight,'user-defined')
         varID = find(w(n_units+1:end));
         for i = varID'
            tw = w(i+n_units);
            if logFlag(i+n_units)
               y = y + tw^2*(exp(x(i))-varNom(i))^2;
               gy(i) = gy(i) + 2*tw^2*(exp(x(i))-varNom(i))*exp(x(i));
            else
               y = y + tw^2*(x(i)-varNom(i))^2;
               gy(i) = gy(i) + 2*tw^2*(x(i)-varNom(i));
            end
         end
      end
   end

   function [y,gy] = min_1N(x)
      y = 0;
      gy = zeros(nVar,1);
      varID = find(w(n_units+1:end));
      for i = varID'
         tw = w(n_units+i);
         if logFlag(i+n_units)
            y = y + tw*exp(x(i));
            gy(i) = tw*exp(x(i));
         else
            y = y + tw*abs(x(i));
            gy(i) = tw*sign(x(i));
         end
      end
   end

   function [c,ceq,g,geq] = neq_F(x)
      c = zeros(2*n_units,1);
      g = zeros(nVar,2*n_units);
      for i = 1:n_units
         l = bds(i,1)-abE(i);
         u = bds(i,2)+abE(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            c(2*i-1,1) = tmpN/tmpD - u;
            c(2*i,1) =  l - tmpN/tmpD;
            grad1 = zeros(nVar,1);
            grad2 = zeros(nVar,1);
            grad1(id1) = (tmpD*N(2:end,2:end)*x(id1)-tmpN*D(2:end,2:end)*x(id1))/tmpD^2;
            grad2(id1) = -(tmpD*N(2:end,2:end)*x(id1)-tmpN*D(2:end,2:end)*x(id1))/tmpD^2;
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id1)]'*quadCoef*[1;x(id1)]-u;
            c(2*i,1) =  l-[1;x(id1)]'*quadCoef*[1;x(id1)];
            grad1 = zeros(nVar,1);
            grad2 = zeros(nVar,1);
            grad1(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            grad2(id1) = -2*quadCoef(2:end,2:end)*x(id1)-2*quadCoef(2:end,1);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      ceq = [];
      geq = [];
   end

end






