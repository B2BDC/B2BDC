function xopt = optParameter_fmincon_withPenalty(obj,id,Hflag,wX,wY,alpha,b2bopt)
% function XOPT = OPTPARAMETER_FMINCON_withPenalty(OBJ,ID,FLAG,HFLAG,WX,WY,ALPHA,B2BOPT) 
% returns the optimal parameter according to the optimization standard

%  Created: Sep 25, 2017    Wenyu Li

if size(wX,1) == 1
   wX = wX';
end
if size(wY,1) == 1
   wY = wY';
end
vList = obj.Variables;
nVar = vList.Length;
nomVal = [vList.Values.NominalValue]';
xbd = vList.calBound;
LB = xbd(:,1);
UB = xbd(:,2);
% if flag && ~Hflag
%    opt = optimoptions('fmincon','Display','none','GradObj','on',...
%       'GradConstr','on','Algorithm','interior-point','Hessian',...
%       'user-supplied','HessFcn',@hessianfcn,'TolFun',1e-6,'TolCon',1e-10);
% else
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','MaxIter',3000,'TolFun',1e-6,'TolCon',1e-10);
% end
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
bds = obj.calBound;
bds = diff(bds,[],2);
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
[xopt,yopt,exitflag] = fmincon(@funxmin,x0,...
   A1,B1,Aeq,Beq,LB,UB,@neq,opt);
if exitflag < 0
   xopt = [];
end

   function [y,gy] = funxmin(x)
      y = 0;
      gy = zeros(nVar,1);
      for i = 1:n_units
         tw = wY(i);
         id1 = idall{i};
         d = units(i).ObservedValue;
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            y = y+tw^2*(tmpN/tmpD-d)^2;
            gy(id1) = gy(id1)+2*tw^2*(tmpN/tmpD-d)/tmpD^2*...
               (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
               tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            y = y + tw^2*(tmpQ - d)^2;
            gy(id1) = gy(id1) + 2*tw^2*(tmpQ-d)*...
               (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
         end
      end
      y = y+alpha*(wX.*x)'*(wX.*x);
      gy = gy + 2*alpha*wX.^2.*x;
   end
   
   function [c,ceq,g,geq] = neq(x)
      if ~Hflag
         c = zeros(2*n_units,1);
         g = zeros(nVar,2*n_units);
         for i = 1:n_units
            l = units(i).LowerBound-abE(i);
            u = units(i).UpperBound+abE(i);
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
      else
         c = [];
         g = [];
      end
      ceq = [];
      geq = [];
   end
   
%    function hessian = hessianfcn(x, lambda)
%       hessian = zeros(nVar);
%       for i = 1:n_units
%          id1 = idall{i};
%          tw = wY(i);
%          quadHessian = 2*Qall{i}(2:end,2:end);
%          hessianMatrix1 = zeros(nVar);
%          hessianMatrix2 = zeros(nVar);
%          hessianMatrix1(id1,id1) = quadHessian;
%          hessianMatrix2(id1,id1) = -quadHessian;
%          hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
%             lambda.ineqnonlin(2*i)*hessianMatrix2 + 2*tw^2*hessianMatrix1;
%       end
%    end
end
