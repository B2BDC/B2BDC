function [yin_result,s,xopt] = absCMfminconNN(obj,disflag,b2bopt)
% subfunction to calculate CM inner bound with fmincon

%  Created: Aug 21, 2018     Wenyu Li

vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
xLB = [vars.LowerBound]';
xUB = [vars.UpperBound]';
xbd = xUB-xLB;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   n3 = size(A0,1);
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
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
bds = obj.calBound;
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
UB = bds(:,2)+abE;
LB = bds(:,1)-abE;
bds = diff(bds,[],2);
allVarnames = obj.VarNames;
idall = cell(n_units,1);
fg = cell(n_units,1);
for j = 1:n_units
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   fg{j} = tmodel.NetGradientFcn;
end
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','MaxIter',5000,...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
if isempty(obj.FeasiblePoint)
   nStart = b2bopt.LocalStart;
   xStart = vList.makeLHSsample(nStart);
   for j = 1:nStart
      [x0,yin_result,exitflag,~,ilam] = fmincon(@funxmin,[0;xStart(j,:)'],...
         A1,B1,Aeq,Beq,[-Inf;xLB],[Inf;xUB],@neq,opt);
      if yin_result<0 && exitflag>0
         break
      end
   end
else
   [x0,yin_result,~,~,ilam] = fmincon(@funxmin,[0;obj.FeasiblePoint],...
      A1,B1,Aeq,Beq,[-Inf;xLB],[Inf;xUB],@neq,opt);
end
xopt = x0(2:end);
yin_result = -yin_result;
s.expu = ilam.ineqnonlin(1:2:2*n_units).*bds;
s.expl = ilam.ineqnonlin(2:2:2*n_units).*bds;
s.varu = ilam.upper(2:end).*xbd;
s.varl = ilam.lower(2:end).*xbd;
if ~isempty(A1)
   s.linu = zeros(n3,1);
   s.linl = zeros(n3,1);
   s.linu(~ieq) = ilam.ineqlin(1:n3).*dA(~ieq);
   s.linl(~ieq) = ilam.ineqlin(n3+1:2*n3).*dA(~ieq);
else
   s.linu = [];
   s.linl = [];
end


   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1;zeros(n_variable,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+1,2*n_units);
      for i = 1:n_units
         id = idall{i}+1;
         c(2*i-1) = units(i).SurrogateModel.eval(x(id)');
         g(id,2*i-1) = fg{i}(x(id));
      end
      c(2:2:end) = -c(1:2:end);
      g(:,2:2:end) = -g(:,1:2:end);
      c(1:2:end) = c(1:2:end)+x(1)-UB;
      c(2:2:end) = c(2:2:end)+LB+x(1);
      g(1,:) = 1;
      ceq = [];
      geq = [];
   end

end

