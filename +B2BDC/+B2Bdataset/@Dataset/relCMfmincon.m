function [yin_result,s,xopt,abE,flag] = relCMfmincon(obj,disflag,b2bopt)
% subfunction to calculate CM inner bound with fmincon

%  Created: Oct 7, 2016     Wenyu Li

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

[idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion;

if obj.ModelDiscrepancyFlag
   MDvar = obj.ModelDiscrepancy.Variables;
   nMD = MDvar.Length;
   MDbd = MDvar.calBound;
   LB = [LB; MDbd(:,1)];
   UB = [UB; MDbd(:,2)];
   if obj.ParameterDiscrepancyFlag
      PDvar = obj.ParameterDiscrepancy.Variables;
      nPD = PDvar.Length;
      PDbd = PDvar.calBound;
      LB = [LB; PDbd(:,1)];
      UB = [UB; PDbd(:,2)];
   else
      nPD = 0;
   end
else
   nMD = 0;
   if obj.ParameterDiscrepancyFlag
      PDvar = obj.ParameterDiscrepancy.Variables;
      nPD = PDvar.Length;
      PDbd = PDvar.calBound;
      LB = [LB; PDbd(:,1)];
      UB = [UB; PDbd(:,2)];
   else
      nPD = 0;
   end
end
xbd = UB-LB;
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
   APD = [zeros(size(APD,1),1) APD];
end

% if flag
%    opt = optimoptions('fmincon','Display','none','GradObj','on',...
%       'GradConstr','on','Algorithm','interior-point','Hessian',...
%       'user-supplied','HessFcn',@hessianfcn,'TolFun',1e-6,'TolCon',1e-6);
% else
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','MaxIter',1500,...
      'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
% end
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
if isempty(obj.FeasiblePoint)
   nStart = b2bopt.LocalStart;
   xStart = zeros(nStart,n_variable+nMD+nPD+1);
   xStart(:,2:n_variable+1) = vList.makeLHSsample(nStart);
   for j = 1:nStart
      [x0,yin_result,~,~,ilam] = fmincon(@funxmin,xStart(j,:)',...
         [A1;APD],[B1;bPD],Aeq,Beq,[-Inf;LB],[Inf;UB],@neq,opt);
      if yin_result<0
         break
      end
   end
else
   xStart = zeros(1,n_variable+nMD+nPD+1);
   xStart(2:n_variable+1) = obj.FeasiblePoint;
   xStart(n_variable+1+1:n_variable+1+nMD) = obj.ModelDiscrepancy.FeasiblePoint;
   xStart(n_variable+1+nMD+1:n_variable+1+nMD+nPD) = obj.ParameterDiscrepancy.FeasiblePoint;
   [x0,yin_result,~,~,ilam] = fmincon(@funxmin,xStart,...
      [A1;APD],[B1;bPD],Aeq,Beq,[0;LB],[Inf;UB],@neq,opt);
end
xopt = x0(2:n_variable+1);
yin_result = -yin_result;
if nMD > 0
   obj.ModelDiscrepancy.FeasiblePoint = x0(n_variable+2:n_variable+nMD+1);
end
if nPD > 0
   obj.ParameterDiscrepancy.FeasiblePoint = x0(n_variable+nMD+2:end);
end
s.expu = ilam.ineqnonlin(1:2:2*n_units).*bds;
s.expl = ilam.ineqnonlin(2:2:2*n_units).*bds;
s.varu = ilam.upper(2:n_variable+nMD+nPD+1).*xbd;
s.varl = ilam.lower(2:n_variable+nMD+nPD+1).*xbd;
n3 = size(A1,1);
if ~isempty(A1)
   s.linu = zeros(n1+n2,1);
   s.linl = zeros(n1+n2,1);
   s.linu(~ieq) = ilam.ineqlin(1:n3).*dA(~ieq);
   s.linl(~ieq) = ilam.ineqlin(n3+1:2*n3).*dA(~ieq);
else
   s.linu = [];
   s.linl = [];
end
if ~isempty(APD)
   n4 = 0.5*size(APD,1);
   dAPD = bPD(1:2:end)+bPD(2:2:end);
   s.linu = [s.linu ; ilam.ineqlin(2*n3+1:2:2*(n3+n4)).*dAPD];
   s.linl = [s.linl ; ilam.ineqlin(2*n3+2:2:2*(n3+n4)).*dAPD];
end



   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1;zeros(n_variable+nMD+nPD,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD+1,2*n_units);
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
            c(2*i-1,1) = [1;x(id+1)]'*quadCoef*[1;x(id+1)]-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-[1;x(id+1)]'*quadCoef*[1;x(id+1)];
%             grad1 = zeros(n_variable+nMD+nPD+1,1);
%             grad2 = zeros(n_variable+nMD+nPD+1,1);
%             grad1(id+1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
%             grad1(1) = u-d;
%             grad2(id+1) = -2*quadCoef(2:end,2:end)*x(id+1)-2*quadCoef(2:end,1);
%             grad2(1) = d-l;
%             g(:,2*i-1) = grad1;
%             g(:,2*i) = grad2;
            g(id+1,2*i-1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
            g(id+1,2*i) = -g(id+1,2*i-1);
            g(1,2*i-1) = u-d;
            g(1,2*i) = d-l;
         end
      end
      ceq = [];
      geq = [];
   end

   function hessian = hessianfcn(x, lambda)
      hessian = zeros(n_variable+1);
      for i = 1:n_units
         id = idall{i};
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable+1);
         hessianMatrix2 = zeros(n_variable+1);
         hessianMatrix1(id+1,id+1) = quadHessian;
         hessianMatrix2(id+1,id+1) = -quadHessian;
         hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
            lambda.ineqnonlin(2*i)*hessianMatrix2;
      end
   end
end